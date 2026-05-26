program LISA_Evolucion_Temporal
    implicit none

    ! --- DEFINICION DE PRECISION ---
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: PI = 3.14159265358979323846_dp

    ! --- PARAMETROS FISICOS ---
    real(dp) :: c = 1.0_dp
    real(dp) :: L = 8.36_dp           ! Brazo LISA
    real(dp) :: w = 0.628_dp          ! Frecuencia angular
    real(dp) :: eps = 1.0_dp          ! Prefactor
    
    ! --- VALOR DE HUBBLE (Solo Físico) ---
    real(dp) :: H_std   = 2.3e-18_dp 
    
    real(dp) :: ZA      = 5.148e15_dp ! Distancia fija (~0.1 Gpc / 5.148e15 s)
    real(dp) :: T_actual

    ! Variables de control
    integer :: i, j
    real(dp) :: alpha_deg, alpha_rad
    ! Array de 101 para cubrir de 0 a 1000s en pasos de 10s
    real(dp) :: tau_results(101)       
    real(dp) :: dummy_int

    open(unit=30, file='lisa_evolucion_T.dat', status='replace')
    
    ! --- ESCRITURA DE CABECERA AUTOMATIZADA ---
    write(30, '(A10)', advance='no') "# Alpha"
    do j = 0, 100
        write(30, '(A6,I4,A1)', advance='no') " T=Z+", j*10, "s"
    end do
    write(30, *) ! Nueva línea tras la cabecera

    print *, "=========================================================="
    print *, " CALCULO DE EVOLUCION TEMPORAL (T > Z) PARA H_STD"
    print *, " Resolucion: 10s | Total casos: 101"
    print *, "=========================================================="

    ! --- BARRIDO ANGULAR (0 a 90 grados) ---
    do i = 0, 90
        alpha_deg = real(i, dp)
        alpha_rad = alpha_deg * (PI / 180.0_dp)

        ! --- BARRIDO TEMPORAL (Pasos de 10s) ---
        do j = 0, 100
            T_actual = ZA + real(j * 10, dp)
            
            call Calcular_Valores_Literales(alpha_rad, H_std, ZA, T_actual, L, w, c, eps, &
                                            dummy_int, tau_results(j+1))
        end do

        ! Escritura de Alpha + 101 columnas de resultados
        write(30, '(F10.2, 101E20.10)') alpha_deg, tau_results
        
        if (mod(i, 10) == 0) print *, "Procesado: ", i, " grados..."
    end do

    close(30)
    print *, "Calculo finalizado. Datos en 'lisa_evolucion_T.dat'"

contains

subroutine Calcular_Valores_Literales(alpha, H, Z, T, L, w, c, eps, v_int, v_tau)
    real(dp), intent(in) :: alpha, H, Z, T, L, w, c, eps
    real(dp), intent(out) :: v_int, v_tau
    
    integer :: j, n_int
    real(dp) :: x, dx, cos_a, Delta_T_val
    real(dp) :: g_val, f_slope, f_const, h_quad, h_lin, h_c
    real(dp) :: val_f, val_g, val_h, R_x, T_x, Theta, amp_func, integral_sum, term_val
    
    n_int = 100000
    dx = 1.0_dp / real(n_int, dp)
    cos_a = cos(alpha)
    
    Delta_T_val = T - Z

    g_val = - (w * H / 2.0_dp)
    f_slope = - (w * H * L)
    f_const = - (w * H * Delta_T_val)
    h_quad = (w * H * (L**2) / 2.0_dp) * (cos_a - 2.0_dp) * cos_a
    h_lin  = L * w * ( 1.0_dp - (H * Delta_T_val + 1.0_dp) * cos_a )
    h_c    = w * Delta_T_val

    integral_sum = 0.0_dp
    do j = 0, n_int
        x = -1.0_dp + real(j, dp) * dx
        
        R_x = Z + x * L * cos_a
        T_x = T + x * L
        
        val_g = g_val
        val_f = f_slope * x + f_const
        val_h = h_quad * x**2 + h_lin * x + h_c
        
        Theta = Z * val_f + (Z**2) * val_g + val_h
        
        amp_func = (1.0_dp + H * T_x) / R_x
        term_val = amp_func * cos(Theta)

        if (j == 0 .or. j == n_int) then
            integral_sum = integral_sum + 0.5_dp * term_val
        else
            integral_sum = integral_sum + term_val
        end if
    end do
    
    v_int = integral_sum * dx
    v_tau = eps * (sin(alpha)**2) * abs(v_int)

end subroutine Calcular_Valores_Literales
end program LISA_Evolucion_Temporal