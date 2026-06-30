program saint_venant_1d
    !-----------------------------------------------------------------------
    ! Solver explicito 1D de las ecuaciones de Saint-Venant (flujo en
    ! lamina libre, regimen no permanente) mediante el esquema
    ! predictor-corrector de MacCormack.
    !-----------------------------------------------------------------------
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: g = 9.81_dp   ! aceleracion de la gravedad [m/s^2]

    integer  :: i, n, nx, nt, unit_out
    real(dp) :: dx, dt, A_izq, A_der, courant, vel_max

    real(dp), allocatable :: A(:), Q(:), U(:)
    real(dp), allocatable :: Anew(:), Qnew(:)

    print *, "=== Solver 1D de Saint-Venant (esquema de MacCormack) ==="

    print *, "Ingrese dx (tamano del paso espacial, en metros):"
    read *, dx

    print *, "Ingrese dt (paso temporal, en segundos):"
    read *, dt

    print *, "Ingrese nx (numero de nodos espaciales, entero >= 3):"
    read *, nx

    print *, "Ingrese nt (numero de pasos temporales, entero >= 1):"
    read *, nt

    print *, "Ingrese el area hidraulica inicial a la izquierda (m^2):"
    read *, A_izq

    print *, "Ingrese el area hidraulica inicial a la derecha (m^2):"
    read *, A_der

    ! --- Validacion de los parametros ingresados ---
    if (nx < 3) then
        print *, "ERROR: nx debe ser mayor o igual a 3."
        stop 1
    end if
    if (nt < 1) then
        print *, "ERROR: nt debe ser mayor o igual a 1."
        stop 1
    end if
    if (dx <= 0.0_dp .or. dt <= 0.0_dp) then
        print *, "ERROR: dx y dt deben ser valores positivos."
        stop 1
    end if
    if (A_izq <= 0.0_dp .or. A_der <= 0.0_dp) then
        print *, "ERROR: las areas hidraulicas iniciales deben ser positivas."
        stop 1
    end if

    allocate(A(nx), Q(nx), U(nx), Anew(nx), Qnew(nx))

    call condiciones_iniciales(A, Q, nx, A_izq, A_der)

    ! --- Verificacion previa de la condicion de estabilidad (CFL) ---
    ! El esquema explicito requiere, de forma aproximada:
    !     dt <= dx / (|u| + c),   con c = sqrt(g*A) (celeridad de la onda)
    vel_max = maxval(abs(Q/A) + sqrt(g*A))
    courant = vel_max*dt/dx
    if (courant > 1.0_dp) then
        print *, "ADVERTENCIA: numero de Courant =", courant, &
                 " (> 1). La simulacion puede volverse inestable."
        print *, "Se recomienda reducir dt o aumentar dx."
    else
        print *, "Numero de Courant =", courant, " -> esquema estable."
    end if

    ! --- Bucle temporal explicito ---
    do n = 1, nt
        call resolver_saint_venant(A, Q, Anew, Qnew, nx, dx, dt, g)
        A = Anew
        Q = Qnew
    end do

    ! --- Velocidad final del flujo ---
    U = Q/A

    ! --- Exportacion de resultados a archivo (para graficar en Python) ---
    open(newunit=unit_out, file="resultados.dat", status="replace", action="write")
    write(unit_out, '(A)') "# x(m)        A(m^2)       Q(m^3/s)     U(m/s)"
    do i = 1, nx
        write(unit_out, '(4(F12.5,2X))') (i-1)*dx, A(i), Q(i), U(i)
    end do
    close(unit_out)

    print *, "Simulacion finalizada. Resultados guardados en resultados.dat"

end program saint_venant_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine condiciones_iniciales(A, Q, nx, A_izq, A_der)
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer,  intent(in)  :: nx
    real(dp), intent(in)  :: A_izq, A_der
    real(dp), intent(out) :: A(nx), Q(nx)

    integer :: i

    ! Discontinuidad escalon en la mitad del dominio (tipo "rotura de
    ! presa"): mitad izquierda con A_izq, mitad derecha con A_der.
    ! Caudal inicial nulo en todo el dominio (el sistema parte del reposo).
    do i = 1, nx
        if (i < nx/2) then
            A(i) = A_izq
        else
            A(i) = A_der
        end if
        Q(i) = 0.0_dp
    end do

end subroutine condiciones_iniciales

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine resolver_saint_venant(A, Q, Anew, Qnew, nx, dx, dt, g)
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer,  intent(in)  :: nx
    real(dp), intent(in)  :: dx, dt, g
    real(dp), intent(in)  :: A(nx), Q(nx)
    real(dp), intent(out) :: Anew(nx), Qnew(nx)

    real(dp) :: Ap(nx), Qp(nx)
    real(dp) :: F1(nx), F2(nx)
    integer  :: i

    ! --- Flujos en forma conservativa, a partir del estado actual ---
    do i = 1, nx
        F1(i) = Q(i)
        F2(i) = Q(i)*Q(i)/A(i) + 0.5_dp*g*A(i)*A(i)
    end do

    ! --- Etapa predictor (diferencias adelantadas en el espacio) ---
    do i = 2, nx-1
        Ap(i) = A(i) - dt/dx*(F1(i+1)-F1(i))
        Qp(i) = Q(i) - dt/dx*(F2(i+1)-F2(i))
    end do

    ! --- Condiciones de frontera (Dirichlet: extremos fijos) ---
    Ap(1) = A(1); Ap(nx) = A(nx)
    Qp(1) = Q(1); Qp(nx) = Q(nx)

    ! --- Flujos recalculados con el estado predicho ---
    do i = 1, nx
        F1(i) = Qp(i)
        F2(i) = Qp(i)*Qp(i)/Ap(i) + 0.5_dp*g*Ap(i)*Ap(i)
    end do

    ! --- Etapa corrector (diferencias atrasadas + promedio temporal) ---
    do i = 2, nx-1
        Anew(i) = 0.5_dp*(A(i)+Ap(i) - dt/dx*(F1(i)-F1(i-1)))
        Qnew(i) = 0.5_dp*(Q(i)+Qp(i) - dt/dx*(F2(i)-F2(i-1)))
    end do

    Anew(1) = A(1); Anew(nx) = A(nx)
    Qnew(1) = Q(1); Qnew(nx) = Q(nx)

end subroutine resolver_saint_venant
