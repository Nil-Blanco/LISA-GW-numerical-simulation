program LISA_Full_Integral
    implicit none

    ! --- DEFINICION DE PRECISION ---
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: PI = 3.14159265358979323846_dp

    ! ==========================================
    !     USER PARAMETERS (MODIFICABLES)
    ! ==========================================
    ! Constantes Fisicas y de la Mision
    real(dp) :: c    = 1.0_dp          ! Velocidad luz (unidades naturales)
    real(dp) :: H    = 2.3e-18_dp      ! Hubble (s^-1)
    real(dp) :: ZA   = 1.0296e17_dp    ! Distancia Fuente (s) [1 Gpc]
    real(dp) :: L    = 8.336_dp        ! Brazo LISA (s) [2.5e6 km]
    real(dp) :: w    = 1.0_dp          ! Frecuencia angular (rad/s)
    real(dp) :: eps  = 1.0_dp          ! Amplitud de la onda (normalizada a 1)
    
    ! Suposicion razonable para el tiempo: T_A aprox Z_A
    real(dp) :: TA

    ! Parametros de simulacion
    integer :: n_angles = 5000      ! Resolucion angular
    integer :: n_int    = 2000        ! Pasos de integracion (trapecio)
    ! ==========================================

    ! Variables internas
    integer :: i, j
    real(dp) :: alpha_deg, alpha_rad, d_alpha
    real(dp) :: cos_a, sin_a
    real(dp) :: A, B, C_phase
    real(dp) :: term_R, term_T, term_ampl, theta
    real(dp) :: x, dx
    real(dp) :: integral_real, integral_imag, magnitude, prefactor
    real(dp) :: max_val, max_angle

    ! Inicializacion
    TA = ZA  ! Asumimos que el tiempo de vuelo es la distancia
    max_val = -1.0_dp
    max_angle = 0.0_dp

    ! Abrir archivo de salida
    open(unit=10, file='datos_integral_completa.dat', status='replace')

    print *, "==============================================="
    print *, " CALCULO DE INTEGRAL COMPLETA (EXACTA)"
    print *, " Segun la formula de la imagen proporcionada"
    print *, "==============================================="

    ! Barrido de angulos de 0 a 90 grados
    d_alpha = 90.0_dp / real(n_angles, dp)

    do i = 1, n_angles
        alpha_deg = real(i, dp) * d_alpha
        alpha_rad = alpha_deg * (PI / 180.0_dp)
        
        cos_a = cos(alpha_rad)
        sin_a = sin(alpha_rad)

        ! --- 1. CALCULO DE COEFICIENTES (Formula exacta de la imagen) ---
        ! A = (w*H*L^2 / 2) * (cos(a) - 2) * cos(a)
        A = (w * H * L**2 / 2.0_dp) * (cos_a - 2.0_dp) * cos_a

        ! B = L*w * ( -H*ZA + (-H*TA + H*ZA - 1)*cos_a + 1 )
        ! Nota: H*ZA y H*TA son casi iguales, se cancelan en el parentesis, 
        ! dejando (-1)*cos_a. La formula se reduce a la aproximacion, 
        ! pero aqui la escribimos completa.
        B = L * w * ( -H*ZA + (-H*TA + H*ZA - 1.0_dp)*cos_a + 1.0_dp )

        ! C = (w/2) * (-2*H*TA*ZA + H*ZA^2 + 2*TA - 2*ZA)
        C_phase = (w / 2.0_dp) * (-2.0_dp*H*TA*ZA + H*ZA**2 + 2.0_dp*TA - 2.0_dp*ZA)


        ! --- 2. INTEGRACION NUMERICA (Regla del Trapecio) ---
        ! x va de -1 a 0
        dx = 1.0_dp / real(n_int, dp)
        integral_real = 0.0_dp
        integral_imag = 0.0_dp

        do j = 0, n_int
            x = -1.0_dp + real(j, dp) * dx
            
            ! Variables dependientes de x (Dentro de la integral)
            ! R(x) = ZA + x*L*cos(alpha)
            term_R = ZA + x * L * cos_a
            
            ! T(x) = TA + x*L (con c=1)
            term_T = TA + x * L / c

            ! Fase Theta(x) = A*x^2 + B*x + C
            theta = A*x**2 + B*x + C_phase

            ! Amplitud completa: (1/R(x)) * (1 + H*T(x))
            ! Asumimos epsilon'_ij = 1 para ver la magnitud escalar
            term_ampl = (1.0_dp + H * term_T) / term_R

            ! Suma para la integral (Euler complejo para sacar magnitud)
            if (j == 0 .or. j == n_int) then
                integral_real = integral_real + 0.5_dp * term_ampl * cos(theta)
                integral_imag = integral_imag + 0.5_dp * term_ampl * sin(theta)
            else
                integral_real = integral_real + term_ampl * cos(theta)
                integral_imag = integral_imag + term_ampl * sin(theta)
            end if
        end do
        
        integral_real = integral_real * dx
        integral_imag = integral_imag * dx
        
        ! Magnitud de la integral (sin prefactores aun)
        magnitude = sqrt(integral_real**2 + integral_imag**2)

        ! --- 3. PREFACTORES EXTERNOS ---
        ! Formula: tau = - (L * eps / 2c) * sin^2(alpha) * Integral
        ! Calculamos el valor absoluto total
        prefactor = (L * eps / (2.0_dp * c)) * (sin_a**2)
        
        magnitude = prefactor * magnitude

        ! --- 4. GUARDAR Y BUSCAR MAXIMO ---
        write(10, *) alpha_deg, magnitude

        if (magnitude > max_val) then
            max_val = magnitude
            max_angle = alpha_deg
        end if

    end do

    close(10)

    print *, "-----------------------------------------------"
    print *, "RESULTADO DE LA OPTIMIZACION (INTEGRAL EXACTA)"
    print *, "Angulo Optimo calculado: ", max_angle, " grados"
    print *, "Valor maximo (tau_GW)  : ", max_val
    print *, "-----------------------------------------------"
    print *, "Datos guardados en 'datos_integral_completa.dat'"

end program LISA_Full_Integral