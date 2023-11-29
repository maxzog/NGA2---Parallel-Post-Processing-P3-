module particle_stats
    implicit none

    public :: particle_type
    type :: particle_type
        integer :: id
        real(4) :: pos(3)
        real(4) :: fld(3)
    end type particle_type

contains

    function norm(v) result(norm_val)
        real(4), intent(in) :: v(3)
        real(4) :: norm_val
        norm_val = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
    end function norm

    function cross_product(a, b) result(cross_val)
        real(4), intent(in) :: a(3), b(3)
        real(4) :: cross_val(3)
        cross_val(1) = a(2)*b(3) - a(3)*b(2)
        cross_val(2) = a(3)*b(1) - a(1)*b(3)
        cross_val(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product

    subroutine uu_lt_cond_r(ps, nb, L, dr, uul, uut, s)
        type(particle_type), intent(in) :: ps(:)
        integer, intent(in) :: nb
        real(4), intent(in) :: L
        real(4), intent(out) :: dr
        real(4), intent(out) :: uul(nb), uut(nb), s
        integer :: i, j, ir, c(nb+1)
        real(4) :: r(3), rll(3), rt1(3), rt2(3)
        real(4) :: rl, rt, sc, s_temp

        dr = norm([L/2.0, L/2.0, L/2.0]) / nb
        s = 0.0
        sc = 0.0
        uul = 0.0
        uut = 0.0
        c = 0

        do i = 1, size(ps)
            do j = i, size(ps)
                call par_perp_u(ps(i), ps(j), L, rll, rt2)
                r = ps(j)%pos - ps(i)%pos
                ir = floor(norm(r) / dr) + 1
                if (ir <= nb .and. ps(i)%id /= ps(j)%id) then
                    rl = dot_product(ps(i)%fld, rll) * dot_product(ps(j)%fld, rll)
                    rt = dot_product(ps(i)%fld, rt2) * dot_product(ps(j)%fld, rt2)
                    uul(ir) = uul(ir) + rl
                    uut(ir) = uut(ir) + rt
                    c(ir) = c(ir) + 1
                end if
                if (ps(i)%id == ps(j)%id) then
                    s_temp = dot_product(ps(i)%fld, ps(j)%fld)
                    s = s + s_temp
                    sc = sc + 1.0
                end if
            end do
        end do

        do i = 1, nb
            if (c(i) == 0) then
                c(i) = 1
            end if
            uul(i) = uul(i) / (c(i))
            uut(i) = uut(i) / (c(i))
        end do
        s = s / (3.0 * sc)
    end subroutine uu_lt_cond_r
    subroutine sf_lt_cond_r(ps, nb, L, dr, uul, uut, s)
        type(particle_type), intent(in) :: ps(:)
        integer, intent(in) :: nb
        real(4), intent(in) :: L
        real(4), intent(out) :: dr
        real(4), intent(out) :: uul(nb), uut(nb), s
        integer :: i, j, ir, c(nb+1)
        real(4) :: r(3), rll(3), rt1(3), rt2(3)
        real(4) :: rl, rt, sc, s_temp

        dr = norm([L/2.0, L/2.0, L/2.0]) / nb
        s = 0.0
        sc = 0.0
        uul = 0.0
        uut = 0.0
        c = 0

        do i = 1, size(ps)
            do j = i, size(ps)
                call par_perp_u(ps(i), ps(j), L, rll, rt2)
                r = ps(j)%pos - ps(i)%pos
                ir = floor(norm(r) / dr) + 1
                if (ir <= nb .and. ps(i)%id /= ps(j)%id) then
                    rl = dot_product(ps(i)%fld-ps(j)%fld, rll) * dot_product(ps(i)%fld-ps(j)%fld, rll)
                    rt = dot_product(ps(i)%fld-ps(j)%fld, rt2) * dot_product(ps(i)%fld-ps(j)%fld, rt2)
                    uul(ir) = uul(ir) + rl
                    uut(ir) = uut(ir) + rt
                    c(ir) = c(ir) + 1
                end if
                if (ps(i)%id == ps(j)%id) then
                    s_temp = dot_product(ps(i)%fld, ps(j)%fld)
                    s = s + s_temp
                    sc = sc + 1.0
                end if
            end do
        end do

        do i = 1, nb
            if (c(i) == 0) then
                c(i) = 1
            end if
            uul(i) = uul(i) / ( c(i))
            uut(i) = uut(i) / ( c(i))
        end do
        s = s / (3.0 * sc)
    end subroutine sf_lt_cond_r

    subroutine rdf_lt_cond_r(ps, nb, L, dr, c)
        type(particle_type), intent(in) :: ps(:)
        integer, intent(in) :: nb
        real(4), intent(in) :: L
        real(4), intent(out) :: dr
        integer :: i, j, ir
        real(4) :: r(3), rll(3), rt1(3), rt2(3)
        real(4) :: rl, rt, sc, s_temp
        real(4) :: N, V, rho, c(nb), pi, Vshell
        dr = norm([L/2.0, L/2.0, L/2.0]) / nb
        c = 0
        V=L**3
        pi=3.14159265
        N=size(ps)*(size(ps)-1)/2
        rho=N/V
        do i = 1, size(ps)
            do j = i, size(ps)
                r = get_minr(ps(i)%pos, ps(j)%pos, L)
                ir = floor(norm(r) / dr) + 1
                if (ir <= nb .and. ps(i)%id /= ps(j)%id) then 
                    c(ir) = c(ir) + 1
                end if
            end do
        end do

        do i = 1, nb
          Vshell=1.33333*pi*((dr*i)**3-(dr*(i-1))**3)
          c(i) = c(i)/(Vshell*rho)
        end do
    end subroutine rdf_lt_cond_r

    function get_minr(x1, x2 , L) result(r)
       real(4), intent(in) :: x1(3), x2(3), L
       real(4) :: r(3)
       r=x2-x1
       r(1)=r(1)-L*NINT(r(1)/L)
       r(2)=r(2)-L*NINT(r(2)/L)
       r(3)=r(3)-L*NINT(r(3)/L)
    end function get_minr

    subroutine par_perp_u(p, q, L, rll, rt2)
        type(particle_type), intent(in) :: p, q
        real(4), intent(in) :: L
        real(4), intent(out) :: rll(3), rt2(3)
        real(4) :: r(3), rt1(3)
        integer :: i

        r = q%pos - p%pos
        do i = 1, 3
            if (r(i) > L/2.0) then
                r(i) = r(i) - L
            end if
            if (r(i) < -L/2.0) then
                r(i) = r(i) + L
            end if
        end do

        rll = r / norm(r)
        rt1 = cross_product(rll, p%fld) / norm(cross_product(rll, p%fld))
        rt2 = cross_product(rt1, rll) / norm(cross_product(rt1, rll))
    end subroutine par_perp_u

end module particle_stats

