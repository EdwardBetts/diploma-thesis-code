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
