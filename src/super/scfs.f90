module scfs 
    
contains 

subroutine count_scfs(nof, nob, norf, norb, nebm, nqnu, qnu_s, qnu_of, qnu_ob, modul, ncf, lid, binom, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: nof, nob, norf, norb, nebm, nqnu
    integer(8), intent(in) :: qnu_s(nqnu), qnu_of(nof, nqnu), qnu_ob(nob, nqnu), modul(nqnu)
    integer(8), intent(in) :: binom(nebm + 1, nob + 1)
    integer(8), intent(out) :: lid(ishft(binom(nebm + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(out) :: ncf
    integer(8) :: qnu_1(nqnu), i, j, cff1, nb1(nob - norb), diml
    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    call omp_set_num_threads(num_th)
    if (disp_std == 1) print *, 'Counting configurations start, number of threads :', omp_get_max_threads()
    diml = ishft(binom(nebm + 1, nob - norb + 1) + 1, nof - norf)
    !$omp parallel shared(nof, nob, norf, norb, nebm, nqnu, qnu_s, qnu_of, qnu_ob, modul, ncf, lid, binom, disp_std, diml), &
    !$omp& private(qnu_1, i, j, cff1, nb1)
    !$omp do
    do i = 0, diml - 1
        if (mod(i + 1, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Counting configurations', &
                i + 1, '*', omp_get_max_threads(), '/', diml
        end if 
        qnu_1 = qnu_s
        cff1 = ibits(i, 0_8, nof - norf)
        do j = norf + 1, nof
            if (ibits(cff1, j - norf - 1, 1_8) == 0) cycle
            qnu_1 = qnu_1 - qnu_of(j, :)
        end do
        nb1 = decode_nb(nob - norb, nebm, ishft(i, -(nof - norf)), binom)
        do j = norb + 1, nob
            qnu_1 = qnu_1 - qnu_ob(j, :) * nb1(j - norb)
        end do
        lid(i + 2) = 0
        call count_scfs_rec(nof, nob, norf, norb, nqnu, qnu_1, qnu_of, qnu_ob, modul, lid(i + 2))
    end do
    !$omp end do
    !$omp end parallel
    lid(1) = 0
    do i = 1, diml
        lid(i + 1) = lid(i + 1) + lid(i)
    end do
    ncf = lid(diml + 1)
    if (disp_std == 1) print *, 'Counting configurations finish, total number :', ncf
end subroutine

recursive subroutine count_scfs_rec(nof, nob, nomf, nomb, nqnu, qnu_1, qnu_of, qnu_ob, modul, ct)
    implicit none
    integer(8), intent(in) :: nof, nob, nomf, nomb, nqnu
    integer(8), intent(in) :: qnu_1(nqnu), qnu_of(nof, nqnu), qnu_ob(nob, nqnu), modul(nqnu)
    integer(8), intent(inout) :: ct
    logical :: flag
    integer(8) :: i, cyc1, qnu_c(nqnu), nomf1

    flag = .true. 
    do i = 1, nqnu 
        cyc1 = modul(i)
        if (cyc1 == 1 .and. qnu_1(i) == 0) cycle
        if (cyc1 > 1 .and. mod(qnu_1(i), cyc1) == 0) cycle
        flag = .false.
        if (cyc1 == 1 .and. qnu_1(i) < 0) return
    end do
    if (flag) then 
        ct = ct + 1 
        return
    end if
    
    if (nomf > 0) then 
        do nomf1 = 1, nomf
            call count_scfs_rec(nof, nob, nomf1 - 1, nomb, nqnu, qnu_1 - qnu_of(nomf1, :), qnu_of, qnu_ob, modul, ct) ! fermion at nomf1
        end do
        call count_scfs_rec(nof, nob, 0_8, nomb, nqnu, qnu_1, qnu_of, qnu_ob, modul, ct) ! no more fermions
    else
        if (nomb == 0) return
        qnu_c = qnu_1
        do 
            call count_scfs_rec(nof, nob, 0_8, nomb - 1, nqnu, qnu_c, qnu_of, qnu_ob, modul, ct) ! one more boson at nomb
            qnu_c = qnu_c - qnu_ob(nomb, :)
            do i = 1, nqnu 
                if (modul(i) == 1 .and. qnu_c(i) < 0) return
            end do
        end do
    end if
end subroutine 

subroutine generate_scfs(nof, nob, norf, norb, nebm, nqnu, qnu_s, qnu_of, qnu_ob, modul, ncf, &
    lid, rid, conff, confb, binom, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: nof, nob, norf, norb, nebm, nqnu, ncf
    integer(8), intent(in) :: qnu_s(nqnu), qnu_of(nof, nqnu), qnu_ob(nob, nqnu), modul(nqnu)
    integer(8), intent(in) :: binom(nebm + 1, nob + 1)
    integer(8), intent(in) :: lid(ishft(binom(nebm + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(out) :: rid(ishft(binom(nebm + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(out) :: conff(ncf), confb(ncf)
    integer(8) :: qnu_1(nqnu), i, j, ct, nbt(nob), cfft, diml
    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    call omp_set_num_threads(num_th)
    if (disp_std == 1) print *, 'Generating configurations start'
    rid = 0
    diml = ishft(binom(nebm + 1, nob - norb + 1) + 1, nof - norf)
    !$omp parallel shared(nof, nob, norf, norb, nebm, nqnu, qnu_s, qnu_of, qnu_ob, modul, ncf, lid, rid, binom), &
    !$omp& shared(conff, confb, disp_std, diml) private(qnu_1, i, j, ct, nbt, cfft)
    !$omp do
    do i = 0, diml - 1
        if (mod(i + 1, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Generating configurations', &
                i + 1, '*', omp_get_max_threads(), '/', diml
        end if 
        if (lid(i + 1) == lid(i + 2)) cycle 
        qnu_1 = qnu_s
        cfft = ishft(ibits(i, 0_8, nof - norf), norf)
        do j = norf + 1, nof
            if (ibits(cfft, j - 1, 1_8) == 0) cycle
            qnu_1 = qnu_1 - qnu_of(j, :)
        end do
        nbt(1 : norb) = 0
        nbt(norb + 1 : nob) = decode_nb(nob - norb, nebm, ishft(i, -(nof - norf)), binom)
        do j = norb + 1, nob
            qnu_1 = qnu_1 - qnu_ob(j, :) * nbt(j)
        end do
        ct = lid(i + 1)
        call generate_scfs_rec(nof, nob, norf, norb, nebm, norf, norb, nqnu, qnu_1, qnu_of, qnu_ob, modul, ncf, &
            ct, cfft, nbt, lid, rid, conff, confb, binom)
    end do
    !$omp end do
    !$omp end parallel
    if (disp_std == 1)  print *, 'Generating configurations finish'
end subroutine

recursive subroutine generate_scfs_rec(nof, nob, norf, norb, nebm, nomf, nomb, nqnu, qnu_1, qnu_of, qnu_ob, modul, ncf, &
    ct, cfft, nbt, lid, rid, conff, confb, binom)
    implicit none
    integer(8), intent(in) :: nof, nob, norf, norb, nebm, nomf, nomb, nqnu, nbt(nob), ncf, cfft
    integer(8), intent(in) :: qnu_1(nqnu), qnu_of(nof, nqnu), qnu_ob(nob, nqnu), modul(nqnu)
    integer(8), intent(in) :: binom(nebm + 1, nob + 1)
    integer(8), intent(in) :: lid(ishft(binom(nebm + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(inout) :: rid(ishft(binom(nebm + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(inout) :: ct, conff(ncf), confb(ncf)
    integer(8) :: i, cyc1, qnu_c(nqnu), nbt_1(nob), li, ri, nomf1
    logical :: flag

    flag = .true. 
    do i = 1, nqnu 
        cyc1 = modul(i)
        if (cyc1 == 1 .and. qnu_1(i) == 0) cycle
        if (cyc1 > 1 .and. mod(qnu_1(i), cyc1) == 0) cycle
        flag = .false.
        if (cyc1 == 1 .and. qnu_1(i) < 0) return
    end do
    if (flag) then 
        ct = ct + 1 
        conff(ct) = cfft 
        confb(ct) = encode_nb(nob, nebm, nbt, binom)
        ri = ishft(encode_nb(norb, nebm, nbt(1 : norb), binom(:, 1 : norb + 1)), norf) + ibits(cfft, 0_8, norf) + 1
        if (rid(ri) > 0) return
        li = ishft(encode_nb(nob - norb, nebm, nbt(norb + 1 : nob), binom(:, 1 : nob - norb + 1)), nof - norf) &
            + ibits(cfft, norf, nof - norf) + 1
        rid(ri) = ct - lid(li)
        return
    end if
    
    if (nomf > 0) then 
        do nomf1 = 1, nomf
            call generate_scfs_rec(nof, nob, norf, norb, nebm, nomf1 - 1, nomb, nqnu, qnu_1 - qnu_of(nomf1, :), &
                qnu_of, qnu_ob, modul, ncf, ct, ibset(cfft, nomf1 - 1), nbt, lid, rid, conff, confb, binom) ! fermion at nomf1
        end do
        call generate_scfs_rec(nof, nob, norf, norb, nebm, 0_8, nomb, nqnu, qnu_1, qnu_of, qnu_ob, modul, ncf, &
            ct, cfft, nbt, lid, rid, conff, confb, binom) ! no more fermions
    else
        if (nomb == 0) return
        qnu_c = qnu_1
        nbt_1 = nbt
        do 
            call generate_scfs_rec(nof, nob, norf, norb, nebm, 0_8, nomb - 1, nqnu, qnu_c, qnu_of, qnu_ob, modul, ncf, &
                ct, cfft, nbt_1, lid, rid, conff, confb, binom)
            nbt_1(nomb) = nbt_1(nomb) + 1
            qnu_c = qnu_c - qnu_ob(nomb, :)
            do i = 1, nqnu 
                if (modul(i) == 1 .and. qnu_c(i) < 0) return
            end do
        end do
    end if
end subroutine

function encode_nb(nob, nebm, nb, binom) result(cfb)
    implicit none 
    integer(8), intent(in) :: nob, nebm, nb(nob), binom(nebm + 1, nob + 1)
    integer(8) :: cfb 
    integer(8) :: pos, neb1, nb1(nob)
    pos = nob
    neb1 = nebm
    cfb = 0
    nb1 = nb
    do 
        if (pos == 0) return
        if (nb1(pos) == 0) then 
            pos = pos - 1
            cycle
        end if
        cfb = cfb + binom(neb1 + 1, pos)
        nb1(pos) = nb1(pos) - 1
        neb1 = neb1 - 1
    end do
end function

function decode_nb(nob, nebm, cfb, binom) result(nb)
    implicit none 
    integer(8), intent(in) :: nob, nebm, cfb, binom(nebm + 1, nob + 1)
    integer(8) :: nb(nob)
    integer(8) :: pos, neb1, cfb1
    pos = nob
    nb = 0
    neb1 = nebm
    cfb1 = cfb
    do 
        if (cfb1 == 0) return 
        if (cfb1 < binom(neb1 + 1, pos)) then 
            pos = pos - 1 
            cycle
        end if
        cfb1 = cfb1 - binom(neb1 + 1, pos)
        nb(pos) = nb(pos) + 1
        neb1 = neb1 - 1
    end do
end function

function search_scf(nof, nob, norf, norb, nebm, lid, rid, cff, nb, binom) result(i)
    implicit none
    integer(8), intent(in) :: nof, nob, norf, norb, nebm
    integer(8), intent(in) :: cff, nb(nob), binom(nebm + 1, nob + 1)
    integer(8), intent(in) :: lid(ishft(binom(nebm + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(in) :: rid(ishft(binom(nebm + 1, norb + 1) + 1, norf))
    integer(8) :: i
    integer(8) :: li, ri 
    li = ishft(encode_nb(nob - norb, nebm, nb(norb + 1 : nob), binom(:, 1 : nob - norb + 1)), nof - norf) &
        + ibits(cff, norf, nof - norf) + 1
    ri = ishft(encode_nb(norb, nebm, nb(1 : norb), binom(:, 1 : norb + 1)), norf) + ibits(cff, 0_8, norf) + 1
    i = lid(li) + rid(ri)
end function

end