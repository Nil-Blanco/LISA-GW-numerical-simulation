program LISA_Simplificado
    implicit none

    ! --- DEFINICION DE PRECISION ---
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: PI = 3.14159265358979323846_dp

    ! ==========================================
    !     PARAMETROS FISICOS (Unidades Naturales c=1)
    ! ==========================================
    ! H valor estandar (s^-1) [Fuente: 34]
    real(dp) :: H_std = 2.3e-18_dp      
    
    ! Distancia ZA = 1 Gpc en segundos [Fuente: 36]
    real(dp) :: ZA    = 1.0296e17_dp    
    
    ! Brazo L = 2.5e6 km en segundos [Fuente: 35]
    real(dp) :: L     = 8.336_dp        
    
    ! Frecuencia w = 1 rad/s (Para ver bien la resonancia)
    real(dp) :: w     = 1.0_dp          

    ! Variables de resultado
    real(dp) :: ang_opt_std, max_tau_std
    real(dp) :: ang_opt_zero, max_tau_zero

    print *, "==================================================="
    print *, " OPTIMIZACION LISA (c=1, Z=T, INTEGRAL NUMERICA)"
    print *, "==================================================="

    ! ------------------------------------------------
    ! CASO 1: Universo con Constante Cosmologica (H_std)
    ! ------------------------------------------------
    call Calcular_Maximo(H_std, ZA, L, w, ang_opt_std, max_tau_std)
    
    print *, " "
    print *, "--- CASO 1: H = 2.3e-18 s^-1 (LambdaCDM) ---"
    print *, "Angulo Optimo (grados): ", ang_opt_std
    print *, "Maximo Tau (relativo) : ", max_tau_std
    print *, "Nota: La teoria (B=0) predice aprox 40.25 grados"

    ! ------------------------------------------------
    ! CASO 2: Universo Plano / Estatico (H=0)
    ! ------------------------------------------------
    call Calcular_Maximo(0.0_dp, ZA, L, w, ang_opt_zero, max_tau_zero)

    print *, " "
    print *, "--- CASO 2: H = 0.0 (Universo Plano) ---"
    print *, "Angulo Optimo (grados): ", ang_opt_zero
    print *, "Maximo Tau (relativo) : ", max_tau_zero
    print *, "Nota: Sin H, domina el seno^2 (deberia ser ~90 grados)"
    print *, "==================================================="

contains

    ! =========================================================
    ! SUBRUTINA: BARRIDO ANGULAR
    ! =========================================================
    subroutine Calcular_Maximo(H_in, Z_in, L_in, w_in, ang_opt, val_max)
        real(dp), intent(in) :: H_in, Z_in, L_in, w_in
        real(dp), intent(out) :: ang_opt, val_max
        
        integer :: i
        real(dp) :: alpha_deg, alpha_rad, tau_val
        
        val_max = -1.0_dp
        ang_opt = 0.0_dp

        ! Barrido de 0 a 90 grados con alta resolucion
        do i = 1, 9000
            alpha_deg = real(i, dp) * 0.01_dp
            alpha_rad = alpha_deg * (PI / 180.0_dp)

            tau_val = Integral_Numerica(alpha_rad, H_in, Z_in, L_in, w_in)

            if (tau_val > val_max) then
                val_max = tau_val
                ang_opt = alpha_deg
            end if
        end do
    end subroutine Calcular_Maximo

    ! =========================================================
    ! FUNCION: INTEGRAL (Regla del Trapecio)
    ! =========================================================
    function Integral_Numerica(alpha, H, Z, L, w) result(tau)
        real(dp), intent(in) :: alpha, H, Z, L, w
        real(dp) :: tau
        
        integer :: j, n_int
        real(dp) :: x, dx
        real(dp) :: A, B, C_val, cos_a
        real(dp) :: R_val, T_val, Theta, Ampl
        real(dp) :: int_real, int_imag, term_val_re, term_val_im

        n_int = 1000
        dx = 1.0_dp / real(n_int, dp)
        cos_a = cos(alpha)

        ! --- COEFICIENTES SIMPLIFICADOS (Z=T, c=1) ---
        ! Basado en Source 20, 21 y 27 del PDF

        ! A = (w*H*L^2 / 2) * (cos(a) - 2)*cos(a)
        A = (w * H * L**2 / 2.0_dp) * (cos_a - 2.0_dp) * cos_a
        
        ! B simplificado: Al ser Z=T, los terminos H*Z se cancelan parcialmente
        ! B = L*w * (1 - H*Z - cos(a))
        B = L * w * (1.0_dp - H * Z - cos_a)
        
        ! C simplificado: Al ser Z=T
        ! C = - w * H * Z^2 / 2
        C_val = - (w * H * Z**2) / 2.0_dp

        ! --- BUCLE DE INTEGRACION ---
        int_real = 0.0_dp
        int_imag = 0.0_dp

        do j = 0, n_int
            x = -1.0_dp + real(j, dp) * dx
            
            ! Trayectorias [Source 11, 16]
            R_val = Z + x * L * cos_a
            T_val = Z + x * L  ! Z=T asumido en el origen

            ! Fase [Source 18]
            Theta = A*x**2 + B*x + C_val

            ! Amplitud [Source 13, 54]
            Ampl = (1.0_dp + H * T_val) / R_val

            ! Integrando complejo: Ampl * e^(i*Theta)
            term_val_re = Ampl * cos(Theta)
            term_val_im = Ampl * sin(Theta)

            ! Trapecio
            if (j == 0 .or. j == n_int) then
                int_real = int_real + 0.5_dp * term_val_re
                int_imag = int_imag + 0.5_dp * term_val_im
            else
                int_real = int_real + term_val_re
                int_imag = int_imag + term_val_im
            end if
        end do
        
        int_real = int_real * dx
        int_imag = int_imag * dx

        ! --- RESULTADO FINAL ---
        ! Tau proporcional a: sin^2(alpha) * Magnitud_Integral
        tau = (sin(alpha)**2) * sqrt(int_real**2 + int_imag**2)

    end function Integral_Numerica

end program LISA_Simplificado