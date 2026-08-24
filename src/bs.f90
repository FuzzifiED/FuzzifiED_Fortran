module bs

use cfs

contains

subroutine generate_trs(no, nor, ncf, lid, rid, conf, perm_o, ph_o, fac_o, id_f, phase, num_th)
    implicit none
    integer(8), intent(in) :: no, nor, ncf
    integer(8), intent(in) :: lid(ibset(0_8, no - nor) + 1), rid(ibset(0_8, nor)), conf(ncf)
    integer(8), intent(in) :: perm_o(no), ph_o(no)
    complex(8), intent(in) :: fac_o(no)
    integer(8), intent(out) :: id_f(ncf)
    complex(8), intent(out) :: phase(ncf)
    integer(8), intent(in) :: num_th 

    integer(8) :: i, o, o1, cf, cf1, cfp, ne, sgn, cf0 
    complex(8) :: phi

    cf0 = 0 
    do o = 0, no - 1
        if (ph_o(o + 1) == 1) cf0 = ibset(cf0, o)
    end do
    call omp_set_num_threads(num_th)
    !$omp parallel shared(no, nor, perm_o, fac_o, ncf, lid, rid, conf, id_f, phase, cf0), &
    !$omp& private(i, o, o1, cf, cf1, cfp, ne, sgn, phi)
    !$omp do 
    do i = 1, ncf
        cf = conf(i)
        cf1 = cf0
        ne = 0
        phi = 1.
        cfp = 0
        do o = no - 1, 0, -1
            o1 = perm_o(o + 1) - 1 
            if (ibits(cf, o, 1_8) == 0) cycle
            phi = phi * fac_o(o + 1)
            cf1 = ieor(cf1, ibset(0_8, o1))
            cfp = ieor(cfp, ibits(cf1, 0_8, o1))
        end do
        sgn = 0
        do o = 0, no - 1
            sgn = ieor(sgn, ibits(cfp, o, 1_8))
        end do
        if (sgn == 1) phi = -phi
        id_f(i) = search_cf(no, nor, lid, rid, cf1)
        phase(i) = phi
    end do
    !$omp end do
    !$omp end parallel
end subroutine

subroutine action_trs(no, nor, &
    ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, &
    ncf_f, dim_f, conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, &
    perm_o, ph_o, fac_o, st_d, st_f, num_th)
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

    integer(8), intent(in) :: perm_o(no), ph_o(no)
    complex(8), intent(in) :: fac_o(no)

    complex(8), intent(in) :: st_d(dim_d)
    complex(8), intent(out) :: st_f(dim_f)
    integer(8), allocatable :: id_f(:)
    complex(8), allocatable :: st_f1(:), phase(:)

    integer(8), intent(in) :: num_th 

    integer(8) :: e, g, g1, i, i1, mult
    complex(8) :: val, fac, fac1

    allocate(id_f(ncf_d))
    allocate(phase(ncf_d))
    call generate_trs(no, nor, ncf_d, lid_f, rid_f, conf_d, perm_o, ph_o, fac_o, id_f, phase, num_th)

    call omp_set_num_threads(num_th)
    st_f = 0
    !$omp parallel shared(no, nor, ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, ncf_f, dim_f), &
    !$omp& shared(conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, perm_o, ph_o, fac_o, st_d, st_f, id_f, phase), &
    !$omp& private(g, g1, e, i, i1, val, fac, fac1, st_f1, mult)
    allocate(st_f1(dim_f))
    st_f1 = 0
    !$omp do
    do g = 1, dim_d
        mult = grsz_d(g)
        do e = 1, mult
            i = grel_d(e, g)
            if (i == -1) cycle
            i1 = id_f(i)
            val = phase(i)
            fac = cffac_d(i)

            g1 = cfgr_f(i1)
            if (g1 == -1) cycle
            fac1 = cffac_f(i1)
            st_f1(g1) = st_f1(g1) + st_d(g) * val * conjg(fac1) * fac
        end do
    end do
    !$omp end do 
    !$omp critical 
    st_f = st_f + st_f1 
    !$omp end critical
    deallocate(st_f1)
    !$omp end parallel
    deallocate(id_f)
    deallocate(phase)
end subroutine

subroutine action_trs_mat(no, nor, &
    ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, &
    ncf_f, dim_f, conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, &
    perm_o, ph_o, fac_o, nst, st_d, st_f, num_th)
    use omp_lib
    implicit none
    integer(8), intent(in) :: no, nor, nst 

    integer(8), intent(in) :: ncf_d, dim_d, szz_d
    integer(8), intent(in) :: conf_d(ncf_d), lid_d(ibset(0_8, no - nor)), rid_d(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_d(ncf_d), grel_d(szz_d, dim_d), grsz_d(dim_d)
    complex(8), intent(in) :: cffac_d(ncf_d)

    integer(8), intent(in) :: ncf_f, dim_f, szz_f
    integer(8), intent(in) :: conf_f(ncf_f), lid_f(ibset(0_8, no - nor) + 1), rid_f(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_f(ncf_f), grel_f(szz_f, dim_f), grsz_f(dim_f)
    complex(8), intent(in) :: cffac_f(ncf_f)

    integer(8), intent(in) :: perm_o(no), ph_o(no)
    complex(8), intent(in) :: fac_o(no)

    complex(8), intent(in) :: st_d(nst, dim_d)
    complex(8), intent(out) :: st_f(nst, dim_f)
    integer(8), allocatable :: id_f(:)
    complex(8), allocatable :: st_f1(:, :), phase(:)

    integer(8), intent(in) :: num_th 

    integer(8) :: e, g, g1, i, i1, mult
    complex(8) :: val, fac, fac1

    allocate(id_f(ncf_d))
    allocate(phase(ncf_d))
    call generate_trs(no, nor, ncf_d, lid_f, rid_f, conf_d, perm_o, ph_o, fac_o, id_f, phase, num_th)

    call omp_set_num_threads(num_th)
    st_f = 0
    !$omp parallel shared(no, nor, nst, ncf_d, dim_d, conf_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, ncf_f), &
    !$omp& shared(dim_f, conf_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, perm_o, ph_o, fac_o, st_d, st_f), &
    !$omp& shared(id_f, phase) private(g, g1, e, i, i1, val, fac, fac1, st_f1, mult)
    allocate(st_f1(nst, dim_f))
    st_f1 = 0
    !$omp do
    do g = 1, dim_d
        mult = grsz_d(g)
        do e = 1, mult
            i = grel_d(e, g)
            if (i == -1) cycle
            i1 = id_f(i)
            val = phase(i)
            fac = cffac_d(i)

            g1 = cfgr_f(i1)
            if (g1 == -1) cycle
            fac1 = cffac_f(i1)
            st_f1(g1, :) = st_f1(g1, :) + st_d(g, :) * val * conjg(fac1) * fac
        end do
    end do
    !$omp end do 
    !$omp critical 
    st_f = st_f + st_f1 
    !$omp end critical
    deallocate(st_f1)
    !$omp end parallel
    deallocate(id_f)
    deallocate(phase)
end subroutine

subroutine generate_bs_cfgr(no, nor, ncf, lid, rid, conf, nqnz, qnz_s, cyc, perm_o, ph_o, fac_o, szz, dim, cfgr, cffac, &
        num_th, disp_std)
    implicit none
    integer(8), intent(in) :: no, nor, ncf
    integer(8), intent(in) :: lid(ibset(0_8, no - nor) + 1), rid(ibset(0_8, nor)), conf(ncf)
    integer(8), intent(in) :: nqnz, cyc(nqnz), szz
    complex(8), intent(in) :: qnz_s(nqnz)
    integer(8), intent(in) :: perm_o(no, nqnz), ph_o(no, nqnz)
    complex(8), intent(in) :: fac_o(no, nqnz)
    integer(8), intent(out) :: dim, cfgr(ncf)
    complex(8), intent(out) :: cffac(ncf)
    integer(8) :: id(szz)
    complex(8) :: fac(szz)
    integer(8), allocatable :: transf(:,:)
    complex(8), allocatable :: phase(:,:)
    integer(8) :: i, j, k, sz
    complex(8) :: wt0
    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    allocate(transf(ncf, nqnz))
    allocate(phase(ncf, nqnz))

    if (disp_std == 1) print *, 'Generating transformations'
    do j = 1, nqnz
        if (abs(qnz_s(j)) < 1.d-6) cycle
        call generate_trs(no, nor, ncf, lid, rid, conf, perm_o(:, j), ph_o(:, j), fac_o(:, j), transf(:, j), phase(:, j), num_th)
    end do
    if (disp_std == 1) print *, 'Generating basis cfgr'
    dim = 0
    cfgr = -1
    do i = 1, ncf
        id(1) = i 
        if (cfgr(id(1)) /= -1) cycle 
        fac(1) = 1.
        sz = 1
        do j = 1, nqnz 
            if (qnz_s(j) == 0) cycle 
            do k = 1, sz * (cyc(j) - 1)
                id(k + sz) = transf(id(k), j)
                fac(k + sz) = fac(k) * phase(id(k), j) * qnz_s(j)
            end do
            sz = sz * cyc(j)
        end do
        wt0 = 0
        do j = 1, sz 
            if (id(1) == id(j)) wt0 = wt0 + fac(j)
        end do 
        if (abs(wt0) < 1.d-6) cycle 
        dim = dim + 1
        cfgr(id(1 : sz)) = dim 
        cffac(id(1 : sz)) = fac(1 : sz) * sqrt(wt0 / sz)
    end do

    deallocate(transf)
    deallocate(phase)

end subroutine

subroutine generate_bs_grel(ncf, szz, dim, cfgr, grel, grsz, disp_std)
    implicit none
    integer(8), intent(in) :: ncf, dim, cfgr(ncf), szz
    integer(8), intent(out) :: grel(szz, dim), grsz(dim)
    integer(8) :: i, g
    integer(8), intent(in) :: disp_std

    if (disp_std == 1) print *, 'Generating basis grel, dimension :', dim
    grel = -1 
    grsz = 0 
    do i = 1, ncf 
        g = cfgr(i)
        if (g == -1) cycle 
        grsz(g) = grsz(g) + 1
        grel(grsz(g), g) = i
    end do
    if (disp_std == 1) print *, 'Generating basis finish'
end subroutine

end module