module op 

use cfs

contains 

function hopping(cfd, cff, no, nc, cstr_tm) result(phase)
    implicit none 
    integer(8), intent(in) :: cfd, no, nc, cstr_tm(2 * nc)
    integer(8), intent(out) :: cff
    integer(8) :: phase
    integer(8) :: i, cfp, cfi, o

    cfi = cfd
    cfp = 0
    phase = 1
    do i = nc * 2 - 1, 1, -2 
        if (cstr_tm(i) == -1) cycle
        o = cstr_tm(i + 1) - 1
        if (cstr_tm(i) == ibits(cfi, o, 1_8)) then 
            phase = 0
            return
        end if
        cfi = ieor(cfi, ibset(0_8, o))
        cfp = ieor(cfp, ibits(cfi, 0_8, o))
    end do 
    cff = cfi
    do i = 0, no - 1
        if (ibits(cfp, i, 1_8) == 1) phase = -phase
    end do
end function

subroutine count_op(no, nor, &
    ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, &
    ncf_f, dim_f, conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, &
    ntm, nc, cstr_tms, fac_tms, red_q, sym_q, nel, colptr, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: no, nor 

    integer(8), intent(in) :: ncf_d, dim_d, szz_d
    integer(8), intent(in) :: conf_d(ncf_d), lid_d(ibset(0_8, no - nor) + 1), rid_d(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_d(ncf_d), grel_d(szz_d, dim_d), grsz_d(dim_d)
    complex(8), intent(in) :: cffac_d(ncf_d)

    integer(8), intent(in) :: ncf_f, dim_f, szz_f
    integer(8), intent(in) :: conf_f(ncf_f), lid_f(ibset(0_8, no - nor) + 1), rid_f(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_f(ncf_f), grel_f(szz_f, dim_f), grsz_f(dim_f)
    complex(8), intent(in) :: cffac_f(ncf_f)

    integer(8), intent(in) :: red_q, sym_q, ntm, nc
    integer(8), intent(in) :: cstr_tms(2 * nc, ntm)
    complex(8), intent(in) :: fac_tms(ntm) ! ComplexF64
    integer(8), intent(out) :: nel, colptr(dim_d + 1)
    
    integer(8) :: g, g1, e, em, i, i1, cf, cf1, t
    integer(8) :: index, index_last, last_id, phase
    integer(8), allocatable :: last_el(:)

    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    call omp_set_num_threads(num_th)
    if (disp_std == 1) print *, 'Counting operator start'
    colptr(1) = 1
    !$omp parallel shared(no, nor, ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, ncf_f, dim_f), &
    !$omp& shared(conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, ntm, nc, cstr_tms, fac_tms, red_q, sym_q, nel), &
    !$omp& shared(colptr, disp_std) private(g, g1, e, em, i, i1, cf, cf1, index, index_last, last_el, last_id, phase, t)
    allocate(last_el(dim_f))
    index = 0
    last_el = 0
    !$omp do
    do g = 1, dim_d
        if (mod(g, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Counting operator, column', g, '*', omp_get_max_threads(), '/', dim_d
        end if 
        index_last = index 
        
        em = 1
        if (red_q == 0) em = szz_d
        do e = 1, em
            i = grel_d(e, g)
            if (i == -1) cycle
            cf = conf_d(i)

            ! 1bd term
            do t = 1, ntm
                phase = hopping(cf, cf1, no, nc, cstr_tms(:, t))
                if (phase == 0) cycle
                i1 = search_cf(no, nor, lid_f, rid_f, cf1)
                if (i1 <= 0 .or. i1 > ncf_f) cycle
                if (cf1 /= conf_f(i1)) cycle

                g1 = cfgr_f(i1)
                if (g1 == -1) cycle 
                if (sym_q > 0 .and. g1 < g) cycle 
                last_id = last_el(g1) 
                if (last_id <= index_last) then 
                    index = index + 1 
                    last_el(g1) = index
                end if
            end do
        end do
        colptr(g + 1) = index - index_last
    end do
    !$omp end do 
    deallocate(last_el)
    !$omp end parallel
    do g = 1, dim_d 
        colptr(g + 1) = colptr(g + 1) + colptr(g)
    end do 
    nel = colptr(dim_d + 1) - 1
    if (disp_std == 1) print *, 'Counting operator finish, element count :', nel
end subroutine

subroutine generate_op(no, nor, &
    ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, &
    ncf_f, dim_f, conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, &
    ntm, nc, cstr_tms, fac_tms, red_q, sym_q, nel, colptr, rowid, elval, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: no, nor 

    integer(8), intent(in) :: ncf_d, dim_d, szz_d
    integer(8), intent(in) :: conf_d(ncf_d), lid_d(ibset(0_8, no - nor)), rid_d(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_d(ncf_d), grel_d(szz_d, dim_d), grsz_d(dim_d)
    complex(8), intent(in) :: cffac_d(ncf_d)

    integer(8), intent(in) :: ncf_f, dim_f, szz_f
    integer(8), intent(in) :: conf_f(ncf_f), lid_f(ibset(0_8, no - nor) + 1), rid_f(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_f(ncf_f), grel_f(szz_f, dim_f), grsz_f(dim_f)
    complex(8), intent(in) :: cffac_f(ncf_f)

    integer(8), intent(in) :: red_q, sym_q, ntm, nc
    integer(8), intent(in) :: cstr_tms(2 * nc, ntm)
    complex(8), intent(in) :: fac_tms(ntm) ! ComplexF64
    integer(8), intent(in) :: nel, colptr(ncf_d + 1)
    integer(8), intent(out) :: rowid(nel)
    complex(8), intent(out) :: elval(nel)

    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    integer(8) :: e, em, g, g1, i, i1, cf, cf1, t
    integer(8) :: index, last_id, phase, mult
    integer(8), allocatable :: last_el(:)
    complex(8) :: val, fac, fac1

    call omp_set_num_threads(num_th)
    if (disp_std == 1) print *, 'Generating operator start'
    !$omp parallel shared(no, nor, ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, ncf_f, dim_f), &
    !$omp& shared(conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, ntm, nc, cstr_tms, fac_tms, red_q, sym_q, nel), &
    !$omp& shared(colptr, rowid, elval, disp_std) private(g, g1, e, em, i, i1, cf, cf1, index, last_el, last_id), &
    !$omp& private(phase, mult, val, fac, fac1, t)
    allocate(last_el(dim_f))
    last_el = 0
    !$omp do
    do g = 1, dim_d
        if (mod(g, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Generating operator column', g, '*', omp_get_max_threads(), '/', dim_d
        end if 
        index = colptr(g) - 1
        mult = grsz_d(g)
        em = 1
        if (red_q == 0) then 
            em = grsz_d(g)
            mult = 1
        end if
        do e = 1, em
            i = grel_d(e, g)
            if (i == -1) cycle
            fac = cffac_d(i)
            cf = conf_d(i)

            do t = 1, ntm
                phase = hopping(cf, cf1, no, nc, cstr_tms(:, t))
                if (phase == 0) cycle
                i1 = search_cf(no, nor, lid_f, rid_f, cf1)
                if (i1 <= 0 .or. i1 > ncf_f) cycle
                if (cf1 /= conf_f(i1)) cycle
                val = fac_tms(t) * phase

                g1 = cfgr_f(i1)
                if (g1 == -1) cycle
                if (sym_q > 0 .and. g1 < g) cycle
                fac1 = cffac_f(i1)
                last_id = last_el(g1)
                if (last_id < colptr(g) .or. last_id >= colptr(g + 1)) then 
                    index = index + 1
                    last_el(g1) = index 
                    rowid(index) = g1
                    elval(index) = val * conjg(fac1) * fac * mult
                else 
                    elval(last_id) = elval(last_id) + val * conjg(fac1) * fac * mult
                end if 
            end do
        end do
    end do
    !$omp end do 
    deallocate(last_el)
    !$omp end parallel
    if (disp_std == 1) print *, 'Generating operator finish'
end subroutine

subroutine generate_op_re(no, nor, &
    ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, &
    ncf_f, dim_f, conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, &
    ntm, nc, cstr_tms, fac_tms, red_q, sym_q, nel, colptr, rowid, elval, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: no, nor 

    integer(8), intent(in) :: ncf_d, dim_d, szz_d
    integer(8), intent(in) :: conf_d(ncf_d), lid_d(ibset(0_8, no - nor)), rid_d(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_d(ncf_d), grel_d(szz_d, dim_d), grsz_d(dim_d)
    complex(8), intent(in) :: cffac_d(ncf_d)

    integer(8), intent(in) :: ncf_f, dim_f, szz_f
    integer(8), intent(in) :: conf_f(ncf_f), lid_f(ibset(0_8, no - nor) + 1), rid_f(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_f(ncf_f), grel_f(szz_f, dim_f), grsz_f(dim_f)
    complex(8), intent(in) :: cffac_f(ncf_f)

    integer(8), intent(in) :: red_q, sym_q, ntm, nc
    integer(8), intent(in) :: cstr_tms(2 * nc, ntm)
    complex(8), intent(in) :: fac_tms(ntm) ! ComplexF64
    integer(8), intent(in) :: nel, colptr(ncf_d + 1)
    integer(8), intent(out) :: rowid(nel)
    real(8), intent(out) :: elval(nel)

    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    integer(8) :: e, em, g, g1, i, i1, cf, cf1, t
    integer(8) :: index, last_id, phase, mult
    integer(8), allocatable :: last_el(:)
    complex(8) :: val, fac, fac1

    call omp_set_num_threads(num_th)
    if (disp_std == 1) print *, 'Generating operator start'
    !$omp parallel shared(no, nor, ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, ncf_f, dim_f), &
    !$omp& shared(conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, ntm, nc, cstr_tms, fac_tms, red_q, sym_q, nel), &
    !$omp& shared(colptr, rowid, elval, disp_std) private(g, g1, e, em, i, i1, cf, cf1, index, last_el, last_id), &
    !$omp& private(phase, mult, val, fac, fac1, t)
    allocate(last_el(dim_f))
    last_el = 0
    !$omp do
    do g = 1, dim_d
        if (mod(g, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Generating operator column', g, '*', omp_get_max_threads(), '/', dim_d
        end if 
        index = colptr(g) - 1
        mult = grsz_d(g)
        em = 1
        if (red_q == 0) then 
            em = grsz_d(g)
            mult = 1
        end if
        do e = 1, em
            i = grel_d(e, g)
            if (i == -1) cycle
            fac = cffac_d(i)
            cf = conf_d(i)

            do t = 1, ntm
                phase = hopping(cf, cf1, no, nc, cstr_tms(:, t))
                if (phase == 0) cycle
                i1 = search_cf(no, nor, lid_f, rid_f, cf1)
                if (i1 <= 0 .or. i1 > ncf_f) cycle
                if (cf1 /= conf_f(i1)) cycle
                val = fac_tms(t) * phase

                g1 = cfgr_f(i1)
                if (g1 == -1) cycle
                if (sym_q > 0 .and. g1 < g) cycle
                fac1 = cffac_f(i1)
                last_id = last_el(g1)
                if (last_id < colptr(g) .or. last_id >= colptr(g + 1)) then 
                    index = index + 1
                    last_el(g1) = index 
                    rowid(index) = g1
                    elval(index) = real(val * conjg(fac1) * fac * mult)
                else 
                    elval(last_id) = elval(last_id) + real(val * conjg(fac1) * fac * mult)
                end if 
            end do
        end do
    end do
    !$omp end do 
    deallocate(last_el)
    !$omp end parallel
    if (disp_std == 1) print *, 'Generating operator finish'
end subroutine

subroutine action_op(no, nor, &
    ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, &
    ncf_f, dim_f, conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, &
    ntm, nc, cstr_tms, fac_tms, red_q, st_d, st_f, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: no, nor 

    integer(8), intent(in) :: ncf_d, dim_d, szz_d
    integer(8), intent(in) :: conf_d(ncf_d), lid_d(ibset(0_8, no - nor)), rid_d(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_d(ncf_d), grel_d(szz_d, dim_d), grsz_d(dim_d)
    complex(8), intent(in) :: cffac_d(ncf_d)

    integer(8), intent(in) :: ncf_f, dim_f, szz_f
    integer(8), intent(in) :: conf_f(ncf_f), lid_f(ibset(0_8, no - nor) + 1), rid_f(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_f(ncf_f), grel_f(szz_f, dim_f), grsz_f(dim_f)
    complex(8), intent(in) :: cffac_f(ncf_f)

    integer(8), intent(in) :: red_q, ntm, nc
    integer(8), intent(in) :: cstr_tms(2 * nc, ntm)
    complex(8), intent(in) :: fac_tms(ntm) ! ComplexF64

    complex(8), intent(in) :: st_d(dim_d)
    complex(8), intent(out) :: st_f(dim_f)
    complex(8), allocatable :: st_f1(:)

    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    integer(8) :: e, em, g, g1, i, i1, cf, cf1, t
    integer(8) :: phase, mult
    complex(8) :: val, fac, fac1

    call omp_set_num_threads(num_th)
    st_f = 0
    if (disp_std == 1) print *, 'Generating operator start'
    !$omp parallel shared(no, nor, ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, ncf_f, dim_f), &
    !$omp& shared(conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, ntm, nc, cstr_tms, fac_tms, red_q, st_d, st_f) &
    !$omp& shared(disp_std) private(g, g1, e, em, i, i1, cf, cf1, phase, val, fac, fac1, t, st_f1, mult)
    allocate(st_f1(dim_f))
    st_f1 = 0
    !$omp do
    do g = 1, dim_d
        if (mod(g, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Generating operator column', g, '*', omp_get_max_threads(), '/', dim_d
        end if 
        mult = grsz_d(g)
        em = 1
        if (red_q == 0) then 
            em = grsz_d(g)
            mult = 1
        end if
        do e = 1, em
            i = grel_d(e, g)
            if (i == -1) cycle
            fac = cffac_d(i)
            cf = conf_d(i)

            do t = 1, ntm
                phase = hopping(cf, cf1, no, nc, cstr_tms(:, t))
                if (phase == 0) cycle
                i1 = search_cf(no, nor, lid_f, rid_f, cf1)
                if (i1 <= 0 .or. i1 > ncf_f) cycle
                if (cf1 /= conf_f(i1)) cycle
                val = fac_tms(t) * phase

                g1 = cfgr_f(i1)
                if (g1 == -1) cycle
                fac1 = cffac_f(i1)
                st_f1(g1) = st_f1(g1) + st_d(g) * val * conjg(fac1) * fac * mult
            end do
        end do
    end do
    !$omp end do 
    !$omp critical 
    st_f = st_f + st_f1 
    !$omp end critical
    deallocate(st_f1)
    !$omp end parallel
    if (disp_std == 1) print *, 'Generating operator finish'
end subroutine

subroutine overlap_op(no, nor, &
    ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, &
    ncf_f, dim_f, conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, &
    ntm, nc, cstr_tms, fac_tms, red_q, st_d, st_f, prod, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: no, nor 

    integer(8), intent(in) :: ncf_d, dim_d, szz_d
    integer(8), intent(in) :: conf_d(ncf_d), lid_d(ibset(0_8, no - nor)), rid_d(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_d(ncf_d), grel_d(szz_d, dim_d), grsz_d(dim_d)
    complex(8), intent(in) :: cffac_d(ncf_d)

    integer(8), intent(in) :: ncf_f, dim_f, szz_f
    integer(8), intent(in) :: conf_f(ncf_f), lid_f(ibset(0_8, no - nor) + 1), rid_f(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_f(ncf_f), grel_f(szz_f, dim_f), grsz_f(dim_f)
    complex(8), intent(in) :: cffac_f(ncf_f)

    integer(8), intent(in) :: red_q, ntm, nc
    integer(8), intent(in) :: cstr_tms(2 * nc, ntm)
    complex(8), intent(in) :: fac_tms(ntm) ! ComplexF64

    complex(8), intent(in) :: st_d(dim_d), st_f(dim_f)
    complex(8), intent(out) :: prod
    complex(8) :: prod1

    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    integer(8) :: e, em, g, g1, i, i1, cf, cf1, t
    integer(8) :: phase, mult
    complex(8) :: val, fac, fac1

    call omp_set_num_threads(num_th)
    prod = 0
    if (disp_std == 1) print *, 'Generating operator start'
    !$omp parallel shared(no, nor, ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, ncf_f, dim_f), &
    !$omp& shared(conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, ntm, nc, cstr_tms, fac_tms, red_q, st_d, st_f), &
    !$omp& shared(prod, disp_std) private(g, g1, e, em, i, i1, cf, cf1, phase, val, fac, fac1, t, prod1, mult)
    prod1 = 0
    !$omp do
    do g = 1, dim_d
        if (mod(g, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Generating operator column', g, '*', omp_get_max_threads()
        end if 
        mult = grsz_d(g)
        em = 1
        if (red_q == 0) then 
            em = grsz_d(g)
            mult = 1
        end if
        do e = 1, em
            i = grel_d(e, g)
            if (i == -1) cycle
            fac = cffac_d(i)
            cf = conf_d(i)

            do t = 1, ntm
                phase = hopping(cf, cf1, no, nc, cstr_tms(:, t))
                if (phase == 0) cycle
                i1 = search_cf(no, nor, lid_f, rid_f, cf1)
                if (i1 <= 0 .or. i1 > ncf_f) cycle
                if (cf1 /= conf_f(i1)) cycle
                val = fac_tms(t) * phase

                g1 = cfgr_f(i1)
                if (g1 == -1) cycle
                fac1 = cffac_f(i1)
                prod1 = prod1 + conjg(st_f(g1)) * st_d(g) * val * conjg(fac1) * fac * mult
            end do
        end do
    end do
    !$omp end do 
    !$omp critical 
    prod = prod + prod1
    !$omp end critical
    !$omp end parallel
    if (disp_std == 1) print *, 'Generating operator finish'
end subroutine

end module