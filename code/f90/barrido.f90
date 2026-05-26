program LISA_Final_Sweep
    implicit none

    ! --- DEFINICION DE PRECISION ---
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: PI = 3.14159265358979323846_dp

    ! --- CONSTANTES FISICAS ---
    real(dp) :: c = 1.0_dp
    real(dp) :: L = 8.336_dp           ! Brazo LISA (s)
    real(dp) :: w = 1.0_dp             ! Frecuencia (rad/s)
    real(dp) :: eps = 1.0_dp           ! Amplitud
    
    ! --- VALORES DE HUBBLE ---
    real(dp) :: H_std  = 2.3e-18_dp    ! H != 0 (LambdaCDM)
    real(dp) :: H_zero = 0.0_dp        ! H = 0 (Estatico)

    ! --- VARIABLES BARRIDO ---
    real(dp) :: Z_Gpc, Z_sec
    real(dp) :: Z_min, Z_max
    real(dp) :: log_min, log_max, log_step, current_log
    real(dp) :: Gpc_to_sec
    integer  :: i, n_steps

    ! --- RESULTADOS ---
    real(dp) :: a0, t0      ! Alpha y Tau para H=0 (Numerico)
    real(dp) :: a1, t1      ! Alpha y Tau para H!=0 (Numerico)
    real(dp) :: a_analit    ! Alpha Analitica (Formula papel)
    real(dp) :: arg_sqrt    ! Argumento raiz cuadrada formula analitica

    ! Factor de conversion: 1 Gpc = 1.0296e17 s
    Gpc_to_sec = 1.0296e17_dp

    ! Configuracion: De 1 kpc (1.0e-6 Gpc) a 100 Gpc
    Z_min = 1.0e-6_dp 
    Z_max = 100.0_dp  
    n_steps = 20      ! Numero de pasos logaritmicos

    log_min = log10(Z_min)
    log_max = log10(Z_max)
    log_step = (log_max - log_min) / real(n_steps - 1, dp)

    print *, "=========================================================================================================="
    print *, "                                   COMPARATIVA LISA: NUMERICA VS ANALITICA                                "
    print *, "=========================================================================================================="

    do i = 1, n_steps
        ! 1. Calcular Z actual (Logaritmico)
        current_log = log_min + real(i-1, dp) * log_step
        Z_Gpc = 10.0_dp**current_log
        Z_sec = Z_Gpc * Gpc_to_sec
        
        ! 2. CALCULO NUMERICO H=0 (Universo Estatico)
        call Encontrar_Maximo(H_zero, Z_sec, Z_sec, L, w, c, eps, a0, t0)

        ! 3. CALCULO NUMERICO H!=0 (Universo LambdaCDM)
        call Encontrar_Maximo(H_std, Z_sec, Z_sec, L, w, c, eps, a1, t1)

        ! 4. CALCULO ANALITICO (Formula del Papel)
        arg_sqrt = (H_std * Z_sec) / 2.0_dp
        
        if (arg_sqrt <= 1.0_dp) then
            a_analit = 2.0_dp * asin(sqrt(arg_sqrt)) * (180.0_dp / PI)
        else
            a_analit = 0.0_dp 
        end if

        ! 5. IMPRIMIR SALIDA
        write(*, '( "Z=", ES9.2, " Gpc -> H=0 => alpha=", F5.2, " \\ Hdif0 => alpha=", F5.2, " \\ Analitica => alpha=", F5.2 )') &
             Z_Gpc, a0, a1, a_analit

    end do
    print *, "=========================================================================================================="

contains

    subroutine Encontrar_Maximo(H_in, Z_in, T_in, L_in, w_in, c_in, eps_in, ang_out, max_out)
        real(dp), intent(in)  :: H_in, Z_in, T_in, L_in, w_in, c_in, eps_in
        real(dp), intent(out) :: ang_out, max_out
        
        integer :: k
        real(dp) :: alpha_deg, alpha_rad, val
        real(dp) :: step
        
        max_out = -1.0_dp
        ang_out = 0.0_dp
        step = 0.1_dp 

        do k = 0, 900
            alpha_deg = real(k, dp) * step
            alpha_rad = alpha_deg * (PI / 180.0_dp)
            val = Integral_Numerica(alpha_rad, H_in, Z_in, T_in, L_in, w_in, c_in, eps_in)
            if (val > max_out) then
                max_out = val
                ang_out = alpha_deg
            end if
        end do
    end subroutine Encontrar_Maximo

    function Integral_Numerica(alpha, H, Z, T, L, w, c, eps) result(val)
        real(dp), intent(in) :: alpha, H, Z, T, L, w, c, eps
        real(dp) :: val
        
        integer :: j, n_int
        real(dp) :: x, dx, cos_a
        real(dp) :: A, B, C_coef
        real(dp) :: R_x, T_x, Theta, Ampl
        real(dp) :: sum_re, sum_im, term_re, term_im
        
        n_int = 500
        dx = 1.0_dp / real(n_int, dp)
        cos_a = cos(alpha)

        A = (w * H * L**2 / 2.0_dp) * (cos_a - 2.0_dp) * cos_a
        B = L * w * ( -H*Z + (-H*T + H*Z - 1.0_dp)*cos_a + 1.0_dp )
        C_coef = (w / 2.0_dp) * (-2.0_dp*H*T*Z + H*Z**2 + 2.0_dp*T - 2.0_dp*Z)

        sum_re = 0.0_dp
        sum_im = 0.0_dp

        do j = 0, n_int
            x = -1.0_dp + real(j, dp) * dx
            R_x = Z + x * L * cos_a
            T_x = T + x * L / c
            
            Theta = A*x**2 + B*x + C_coef
            Ampl = (1.0_dp + H * T_x) / R_x
            
            term_re = Ampl * cos(Theta)
            term_im = Ampl * sin(Theta)

            if (j==0 .or. j==n_int) then
                sum_re = sum_re + 0.5_dp * term_re
                sum_im = sum_im + 0.5_dp * term_im
            else
                sum_re = sum_re + term_re
                sum_im = sum_im + term_im
            end if
        end do
        
        sum_re = sum_re * dx
        sum_im = sum_im * dx
        
        ! Esta linea es la que daba problemas antes por las etiquetas
        val = abs( (L * eps / (2.0_dp * c)) * (sin(alpha)**2) * sqrt(sum_re**2 + sum_im**2) )

    end function Integral_Numerica

end program LISA_Final_Sweep