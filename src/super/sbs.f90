module sbs

use scfs

contains

subroutine generate_strs(nof, nob, norf, norb, nebm, ncf, lid, rid, conff, confb, perm_of, perm_ob, & 
    ph_of, fac_of, fac_ob, id_f, phase, binom, num_th)
    implicit none
    integer(8), intent(in) :: nof, nob, norf, norb, nebm, ncf
    integer(8), intent(in) :: binom(nebm + 1, nob + 1)
    integer(8), intent(in) :: lid(ishft(binom(nebm + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(in) :: rid(ishft(binom(nebm + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(in) :: conff(ncf), confb(ncf)
    integer(8), intent(in) :: perm_of(nof), perm_ob(nob), ph_of(nof)
    complex(8), intent(in) :: fac_of(nof), fac_ob(nob)
    integer(8), intent(out) :: id_f(ncf)
    complex(8), intent(out) :: phase(ncf)
    integer(8), intent(in) :: num_th 

    integer(8) :: i, of, of1, cff, cff1, cffp, sgn, cff0 
    integer(8) :: ob, ob1, cfb, nb(nob), nb1(nob)
    complex(8) :: phi

    cff0 = 0 
    do of = 0, nof - 1
        if (ph_of(of + 1) == 1) cff0 = ibset(cff0, of)
    end do
    call omp_set_num_threads(num_th)
    !$omp parallel shared(nof, nob, ncf, lid, rid, conff, confb, perm_of, perm_ob, ph_of, fac_of, fac_ob, id_f, phase, cff0), &
    !$omp& private(i, of, of1, cff, cff1, cffp, sgn, ob, ob1, cfb, nb, nb1, phi)
    !$omp do 
    do i = 1, ncf
        phi = 1.
        cfb = confb(i)
        nb = decode_nb(nob, nebm, cfb, binom)
        nb1 = 0 
        do ob = 1, nob 
            ob1 = perm_ob(ob)
            phi = phi * fac_ob(ob) ** nb(ob)
            nb1(ob1) = nb(ob)
        end do

        cff = conff(i)
        cff1 = cff0
        cffp = 0
        do of = nof - 1, 0, -1
            of1 = perm_of(of + 1) - 1 
            if (ibits(cff, of, 1_8) == 0) cycle
            phi = phi * fac_of(of + 1)
            cff1 = ieor(cff1, ibset(0_8, of1))
            cffp = ieor(cffp, ibits(cff1, 0_8, of1))
        end do
        sgn = 0
        do of = 0, nof - 1
            sgn = ieor(sgn, ibits(cffp, of, 1_8))
        end do
        if (sgn == 1) phi = -phi
        id_f(i) = search_scf(nof, nob, norf, norb, nebm, lid, rid, cff1, nb1, binom)
        phase(i) = phi
    end do
    !$omp end do
    !$omp end parallel
end subroutine

subroutine action_strs(nof, nob, norf, norb, &
    nebm_d, ncf_d, dim_d, conff_d, confb_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, binom_d, &
    nebm_f, ncf_f, dim_f, conff_f, confb_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, binom_f, &
    perm_of, perm_ob, ph_of, fac_of, fac_ob, st_d, st_f, num_th)
    use omp_lib
    implicit none
    integer(8), intent(in) :: nof, nob, norf, norb

    integer(8), intent(in) :: nebm_d, ncf_d, dim_d, szz_d
    integer(8), intent(in) :: conff_d(ncf_d), confb_d(ncf_d)
    integer(8), intent(in) :: binom_d(nebm_d + 1, nob + 1)
    integer(8), intent(in) :: lid_d(ishft(binom_d(nebm_d + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(in) :: rid_d(ishft(binom_d(nebm_d + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(in) :: cfgr_d(ncf_d), grel_d(szz_d, dim_d), grsz_d(dim_d)
    complex(8), intent(in) :: cffac_d(ncf_d)

    integer(8), intent(in) :: nebm_f, ncf_f, dim_f, szz_f
    integer(8), intent(in) :: conff_f(ncf_f), confb_f(ncf_f)
    integer(8), intent(in) :: binom_f(nebm_f + 1, nob + 1)
    integer(8), intent(in) :: lid_f(ishft(binom_f(nebm_f + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(in) :: rid_f(ishft(binom_f(nebm_f + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(in) :: cfgr_f(ncf_f), grel_f(szz_f, dim_f), grsz_f(dim_f)
    complex(8), intent(in) :: cffac_f(ncf_f)

    integer(8), intent(in) :: perm_of(nof), perm_ob(nob), ph_of(nof)
    complex(8), intent(in) :: fac_of(nof), fac_ob(nob)
    
    complex(8), intent(in) :: st_d(dim_d)
    complex(8), intent(out) :: st_f(dim_f)
    integer(8), intent(in) :: num_th 
    integer(8), allocatable :: id_f(:)
    complex(8), allocatable :: st_f1(:), phase(:)

    integer(8) :: e, g, g1, i, i1, mult
    complex(8) :: val, fac, fac1

    allocate(id_f(ncf_d))
    allocate(phase(ncf_d))
    call generate_strs(nof, nob, norf, norb, nebm_f, ncf_d, lid_f, rid_f, conff_d, confb_d, & 
        perm_of, perm_ob, ph_of, fac_of, fac_ob, id_f, phase, binom_f, num_th)

    call omp_set_num_threads(num_th)
    st_f = 0
    !$omp parallel shared(nof, nob, norf, norb), &
    !$omp& shared(nebm_d, ncf_d, dim_d, conff_d, confb_d, lid_d, rid_d, szz_d, cfgr_d, cffac_d, grel_d, grsz_d, binom_d), &
    !$omp& shared(nebm_f, ncf_f, dim_f, conff_f, confb_f, lid_f, rid_f, szz_f, cfgr_f, cffac_f, grel_f, grsz_f, binom_f), &
    !$omp& shared(perm_of, perm_ob, ph_of, fac_of, fac_ob, st_d, st_f, id_f, phase), & 
    !$omp& private(g, g1, e, i, i1, mult, val, fac, fac1, st_f1)
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

subroutine generate_sbs_cfgr(nof, nob, norf, norb, nebm, ncf, lid, rid, conff, confb, nqnz, qnz_s, cyc, &
        perm_of, perm_ob, ph_of, fac_of, fac_ob, szz, dim, cfgr, cffac, binom, num_th, disp_std)
    implicit none
    integer(8), intent(in) :: nof, nob, norf, norb, nebm, ncf
    integer(8), intent(in) :: binom(nebm + 1, nob + 1)
    integer(8), intent(in) :: lid(ishft(binom(nebm + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(in) :: rid(ishft(binom(nebm + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(in) :: conff(ncf), confb(ncf)
    integer(8), intent(in) :: nqnz, cyc(nqnz), szz
    complex(8), intent(in) :: qnz_s(nqnz)
    integer(8), intent(in) :: perm_of(nof, nqnz), perm_ob(nob, nqnz), ph_of(nof, nqnz)
    complex(8), intent(in) :: fac_of(nof, nqnz), fac_ob(nob, nqnz)
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
        call generate_strs(nof, nob, norf, norb, nebm, ncf, lid, rid, conff, confb, perm_of(:, j), perm_ob(:, j), ph_of(:, j), &
            fac_of(:, j), fac_ob(:, j), transf(:, j), phase(:, j), binom, num_th)
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

subroutine generate_sbs_grel(ncf, szz, dim, cfgr, grel, grsz, disp_std)
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