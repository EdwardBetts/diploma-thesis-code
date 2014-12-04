subroutine fsum(inarr, dim, vout)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(out) :: vout
    real(8), intent(in), dimension(dim) :: inarr

    vout = sum(inarr)
end subroutine

subroutine fmin(inarr, dim, vout)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(out) :: vout
    real(8), intent(in), dimension(dim) :: inarr

    vout = minval(inarr)
end subroutine

subroutine trigger(thres, inarr, range, leninarr, vout)
    implicit none
    integer :: i
    integer, intent(in) :: thres, leninarr
    integer, intent(in), dimension(2) :: range
    integer, intent(in), dimension(leninarr) :: inarr
    integer, intent(out), dimension(200) :: vout

    do i = range(1) - 1, range(2) - 1
        if ((inarr(i) < thres) .and. (inarr(i + 20) < thres)) then
            vout = inarr(i - 20:i + 180)
            exit
        end if
    end do
end subroutine

subroutine trig_ind_pos(thres, inarr, start, dim, vout)
    implicit none
    integer :: i
    integer, intent(in) :: dim, start
    real(8), intent(in) :: thres
    real(8), intent(in), dimension(dim) :: inarr
    integer, intent(out) :: vout

    vout = start
    do i = start, dim - 6
        if ((inarr(i + 1) > thres) .and. (inarr(i + 6) > thres) .and. (inarr(i - 4) < thres)) then
            vout = i
            exit
        end if
    end do
end subroutine

subroutine trig_ind_neg(thres, inarr, start, dim, vout)
    implicit none
    integer :: i
    integer, intent(in) :: dim, start
    real(8), intent(in) :: thres
    real(8), intent(in), dimension(dim) :: inarr
    integer, intent(out) :: vout

    vout = start
    do i = start, dim - 6
        if ((inarr(i + 1) < thres) .and. (inarr(i + 6) < thres) .and. (inarr(i - 4) > thres)) then
            vout = i
            exit
        end if
    end do
end subroutine

subroutine max_and_ind(inarr, dim, start, vout)
    implicit none
    integer :: i
    integer, intent(in) :: dim, start
    real(8) :: maxv
    real(8), intent(in), dimension(dim) :: inarr
    real(8), intent(out), dimension(2) :: vout

    vout(2) = start
    maxv = maxval(inarr(start + 1:start + 31))
    vout(1) = maxv
    do i = start, start + 30
        if (inarr(i + 1) .eq. maxv) then
            vout(2) = real(i)
            exit
        end if
    end do
end subroutine

subroutine min_and_ind(inarr, dim, start, vout)
    implicit none
    integer :: i
    integer, intent(in) :: dim, start
    real(8) :: minv
    real(8), intent(in), dimension(dim) :: inarr
    real(8), intent(out), dimension(2) :: vout

    vout(2) = start
    minv = minval(inarr(start + 1:start + 31))
    vout(1) = minv
    do i = start, start + 30
        if (inarr(i + 1) .eq. minv) then
            vout(2) = real(i)
            exit
        end if
    end do
end subroutine
