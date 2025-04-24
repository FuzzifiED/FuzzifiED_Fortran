module cfs_neg
    
contains

subroutine count_cfs_neg(no, nor, nqnu, qnu_s, qnu_o, modul, ncf, lid, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: no, nor, nqnu
    integer(8), intent(in) :: qnu_s(nqnu), qnu_o(no, nqnu), modul(nqnu)
    integer(8), intent(out) :: lid(ibset(0_8, no - nor) + 1)
    integer(8), intent(out) :: ncf
    integer(8) :: qnu_1(nqnu), i, j
    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    call omp_set_num_threads(num_th)
    if (disp_std == 1) print *, 'Counting configurations start, number of threads :', omp_get_max_threads()
    !$omp parallel shared(no, nor, nqnu, qnu_s, qnu_o, modul, ncf, lid, disp_std) private(qnu_1, i, j)
    !$omp do
    do i = 0, ibset(0_8, no - nor) - 1
        if (mod(i + 1, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Counting configurations', &
                i + 1, '*', omp_get_max_threads(), '/', ibset(0_8, no - nor)
        end if 
        qnu_1 = qnu_s
        do j = nor + 1, no
            if (ibits(i, j - nor - 1, 1_8) == 0) cycle
            qnu_1 = qnu_1 - qnu_o(j, :)
        end do
        lid(i + 2) = 0
        call count_cfs_neg_rec(no, nor, nqnu, qnu_1, qnu_o, modul, lid(i + 2))
    end do
    !$omp end do
    !$omp end parallel
    lid(1) = 0
    do i = 1, ibset(0_8, no - nor)
        lid(i + 1) = lid(i + 1) + lid(i)
    end do
    ncf = lid(ibset(0_8, no - nor) + 1)
    if (disp_std == 1) print *, 'Counting configurations finish, total number :', ncf
end subroutine

recursive subroutine count_cfs_neg_rec(no, nom, nqnu, qnu_1, qnu_o, modul, ct)
    implicit none
    integer(8), intent(in) :: no, nom, nqnu
    integer(8), intent(in) :: qnu_1(nqnu), qnu_o(no, nqnu), modul(nqnu)
    integer(8), intent(inout) :: ct
    logical :: flag
    integer(8) :: i, nom1, cyc1

    flag = .true. 
    do i = 1, nqnu
        cyc1 = modul(i)
        if (cyc1 == 1 .and. qnu_1(i) == 0) cycle
        if (cyc1 > 1 .and. mod(qnu_1(i), cyc1) == 0) cycle
        flag = .false.
    end do
    if (flag) ct = ct + 1 
    do nom1 = 1, nom
        call count_cfs_neg_rec(no, nom1 - 1, nqnu, qnu_1 - qnu_o(nom1, :), qnu_o, modul, ct)
    end do
end subroutine

subroutine generate_cfs_neg(no, nor, nqnu, qnu_s, qnu_o, modul, ncf, lid, rid, conf, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: no, nor, nqnu, ncf
    integer(8), intent(in) :: qnu_s(nqnu), qnu_o(no, nqnu), modul(nqnu)
    integer(8), intent(in) :: lid(ibset(0_8, no - nor) + 1)
    integer(8), intent(out) :: conf(ncf), rid(ibset(0_8, nor))
    integer(8) :: qnu_1(nqnu), i, j, ct
    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    call omp_set_num_threads(num_th)
    if (disp_std == 1) print *, 'Generating configurations start'
    rid = 0
    !$omp parallel shared(no, nor, nqnu, qnu_s, qnu_o, modul, ncf, lid, rid, conf, disp_std) private(qnu_1, i, j, ct)
    !$omp do
    do i = 0, ibset(0_8, no - nor) - 1
        if (mod(i + 1, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Generating configurations', &
                i + 1, '*', omp_get_max_threads(), '/', ibset(0_8, no - nor)
        end if 
        if (lid(i + 1) == lid(i + 2)) cycle 
        qnu_1 = qnu_s
        do j = nor + 1, no
            if (ibits(i, j - nor - 1, 1_8) == 0) cycle
            qnu_1 = qnu_1 - qnu_o(j, :)
        end do
        ct = lid(i + 1)
        call generate_cfs_neg_rec(no, nor, nor, nqnu, qnu_1, qnu_o, modul, ncf, ct, ishft(i, nor), lid, rid, conf)
    end do
    !$omp end do
    !$omp end parallel
    if (disp_std == 1)  print *, 'Generating configurations finish'
end subroutine

recursive subroutine generate_cfs_neg_rec(no, nor, nom, nqnu, qnu_1, qnu_o, modul, ncf, ct, tmp, lid, rid, conf)
    implicit none
    integer(8), intent(in) :: no, nor, nom, nqnu, tmp, ncf
    integer(8), intent(in) :: qnu_1(nqnu), qnu_o(no, nqnu), modul(nqnu)
    integer(8), intent(in) :: lid(ibset(0_8, no - nor) + 1)
    integer(8), intent(inout) :: ct, conf(ncf), rid(ibset(0_8, nor))
    logical :: flag
    integer(8) :: i, nom1, li, ri, cyc1

    flag = .true. 
    do i = 1, nqnu
        cyc1 = modul(i)
        if (cyc1 == 1 .and. qnu_1(i) == 0) cycle
        if (cyc1 > 1 .and. mod(qnu_1(i), cyc1) == 0) cycle
        flag = .false.
    end do
    if (flag) then 
        ct = ct + 1 
        conf(ct) = tmp 
        ri = ibits(tmp, 0_8, nor)
        if (rid(ri + 1) == 0) then
            li = ibits(tmp, nor, no - nor)
            rid(ri + 1) = ct - lid(li + 1)
        end if
    end if
    do nom1 = 1, nom
        call generate_cfs_neg_rec(no, nor, nom1 - 1, nqnu, qnu_1 - qnu_o(nom1, :), qnu_o, &
            modul, ncf, ct, ibset(tmp, nom1 - 1), lid, rid, conf)
    end do
end subroutine

end module
