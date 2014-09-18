subroutine fsum(inarr, dim, vout)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(out) :: vout
    real(8), intent(in), dimension(dim) :: inarr

    vout = sum(inarr)
end subroutine
