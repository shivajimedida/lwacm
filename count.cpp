
void alpha_0_call()
{      
      p_load = p[t_now][x][y][z];
      
      u_load[0] = u[t_now][x][y][z][0];
      u_load[1] = u[t_now][x][y][z][1];
      u_load[2] = u[t_now][x][y][z][2];
      
      u_xi = 0;
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = 0;
      
      f_e = 0.33333 * p_load * ( 1 - 1.5*u_2 );
      
      f_e_o =  0;
      
      f_x_t =  0;
      
      f[0] =  f_e;
}

void alpha_1_call()
{
      p_load = p[t_now][x-1][y][z];
      
      u_load[0] = u[t_now][x-1][y][z][0];
      u_load[1] = u[t_now][x-1][y][z][1];
      u_load[2] = u[t_now][x-1][y][z][2];
      
      u_xi = u_load[0];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0];
      
      
      f_e = 0.055556 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.16667 * p_load * u_xi;
      
      f_x_t = 0.16667 * p[t_now][x][y][z] * u_x_xi;
      
      f[1] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_2_call()
{
      p_load = p[t_now][x+1][y][z];
      
      u_load[0] = u[t_now][x+1][y][z][0];
      u_load[1] = u[t_now][x+1][y][z][1];
      u_load[2] = u[t_now][x+1][y][z][2];
      
      u_xi = u_load[0]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] * (-1);
      
         
      f_e = 0.055556 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.16667 * p_load * u_xi;
      
      f_x_t = 0.16667 * p[t_now][x][y][z] * u_x_xi;
      
      f[2] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_3_call()
{
      p_load = p[t_now][x][y-1][z];
      
      u_load[0] = u[t_now][x][y-1][z][0];
      u_load[1] = u[t_now][x][y-1][z][1];
      u_load[2] = u[t_now][x][y-1][z][2];
      
      u_xi = u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1];
      
      
      f_e = 0.055556 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.16667 * p_load * u_xi;
      
      f_x_t = 0.16667 * p[t_now][x][y][z] * u_x_xi;
      
      f[3] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_4_call()
{
      p_load = p[t_now][x][y+1][z];
      
      u_load[0] = u[t_now][x][y+1][z][0];
      u_load[1] = u[t_now][x][y+1][z][1];
      u_load[2] = u[t_now][x][y+1][z][2];
      
      u_xi = u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1] * (-1);
      
         
      f_e = 0.055556 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.16667 * p_load * u_xi;
      
      f_x_t = 0.16667 * p[t_now][x][y][z] * u_x_xi;
      
      f[4] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_5_call()
{
      p_load = p[t_now][x][y][z-1];
      
      u_load[0] = u[t_now][x][y][z-1][0];
      u_load[1] = u[t_now][x][y][z-1][1];
      u_load[2] = u[t_now][x][y][z-1][2];
      
      u_xi = u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][2];
      
      
      f_e = 0.055556 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.16667 * p_load * u_xi;
      
      f_x_t = 0.16667 * p[t_now][x][y][z] * u_x_xi;
      
      f[5] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_6_call()
{
      p_load = p[t_now][x][y][z+1];
      
      u_load[0] = u[t_now][x][y][z+1][0];
      u_load[1] = u[t_now][x][y][z+1][1];
      u_load[2] = u[t_now][x][y][z+1][2];
      
      u_xi = u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][2] * (-1);
      
         
      f_e = 0.055556 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.16667 * p_load * u_xi;
      
      f_x_t = 0.16667 * p[t_now][x][y][z] * u_x_xi;
      
      f[6] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_7_call()
{
      p_load = p[t_now][x-1][y-1][z];
      
      u_load[0] = u[t_now][x-1][y-1][z][0];
      u_load[1] = u[t_now][x-1][y-1][z][1];
      u_load[2] = u[t_now][x-1][y-1][z][2];
      
      u_xi = u_load[0] + u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] + u[t_now][x][y][z][1];
      
      
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[7] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_8_call()
{
      p_load = p[t_now][x+1][y-1][z];
      
      u_load[0] = u[t_now][x+1][y-1][z][0];
      u_load[1] = u[t_now][x+1][y-1][z][1];
      u_load[2] = u[t_now][x+1][y-1][z][2];
      
      u_xi = u_load[0]*(-1) + u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] * (-1) + u[t_now][x][y][z][1];
      
         
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[8] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_9_call()
{
      p_load = p[t_now][x-1][y+1][z];
      
      u_load[0] = u[t_now][x-1][y+1][z][0];
      u_load[1] = u[t_now][x-1][y+1][z][1];
      u_load[2] = u[t_now][x-1][y+1][z][2];
      
      u_xi = u_load[0] + u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] + u[t_now][x][y][z][1] * (-1);
      
      
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[9] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_10_call()
{
      p_load = p[t_now][x+1][y+1][z];
      
      u_load[0] = u[t_now][x+1][y+1][z][0];
      u_load[1] = u[t_now][x+1][y+1][z][1];
      u_load[2] = u[t_now][x+1][y+1][z][2];
      
      u_xi = u_load[0]*(-1) + u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] * (-1) + u[t_now][x][y][z][1] * (-1);
      
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[10] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_11_call()
{
      p_load = p[t_now][x-1][y][z-1];
      
      u_load[0] = u[t_now][x-1][y][z-1][0];
      u_load[1] = u[t_now][x-1][y][z-1][1];
      u_load[2] = u[t_now][x-1][y][z-1][2];
      
      u_xi = u_load[0] + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] + u[t_now][x][y][z][2];
      
      
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[11] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_12_call()
{
      p_load = p[t_now][x+1][y][z-1];
      
      u_load[0] = u[t_now][x+1][y][z-1][0];
      u_load[1] = u[t_now][x+1][y][z-1][1];
      u_load[2] = u[t_now][x+1][y][z-1][2];
      
      u_xi = u_load[0]*(-1) + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] * (-1) + u[t_now][x][y][z][2];
      
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[12] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_13_call()
{
      p_load = p[t_now][x-1][y][z+1];
      
      u_load[0] = u[t_now][x-1][y][z+1][0];
      u_load[1] = u[t_now][x-1][y][z+1][1];
      u_load[2] = u[t_now][x-1][y][z+1][2];
      
      u_xi = u_load[0] + u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] + u[t_now][x][y][z][2] * (-1);
      
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[13] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_14_call()
{
      p_load = p[t_now][x+1][y][z+1];
      
      u_load[0] = u[t_now][x+1][y][z+1][0];
      u_load[1] = u[t_now][x+1][y][z+1][1];
      u_load[2] = u[t_now][x+1][y][z+1][2];
      
      u_xi = u_load[0]*(-1) + u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] * (-1) + u[t_now][x][y][z][2] * (-1);
       
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[14] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_15_call()
{
      p_load = p[t_now][x][y-1][z-1];
      
      u_load[0] = u[t_now][x][y-1][z-1][0];
      u_load[1] = u[t_now][x][y-1][z-1][1];
      u_load[2] = u[t_now][x][y-1][z-1][2];
      
      u_xi = u_load[1] + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1] + u[t_now][x][y][z][2];
      
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[15] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_16_call()
{
      p_load = p[t_now][x][y+1][z-1];
      
      u_load[0] = u[t_now][x][y+1][z-1][0];
      u_load[1] = u[t_now][x][y+1][z-1][1];
      u_load[2] = u[t_now][x][y+1][z-1][2];
      
      u_xi = u_load[1]*(-1) + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1] * (-1) + u[t_now][x][y][z][2];
      
         
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[16] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_17_call()
{
      p_load = p[t_now][x][y-1][z+1];
      
      u_load[0] = u[t_now][x][y-1][z+1][0];
      u_load[1] = u[t_now][x][y-1][z+1][1];
      u_load[2] = u[t_now][x][y-1][z+1][2];
      
      u_xi = u_load[1] + u_load[2] * (-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1] + u[t_now][x][y][z][2] * (-1);
      
      
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[17] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}

void alpha_18_call()
{
      p_load = p[t_now][x][y+1][z+1];
      
      u_load[0] = u[t_now][x][y+1][z+1][0];
      u_load[1] = u[t_now][x][y+1][z+1][1];
      u_load[2] = u[t_now][x][y+1][z+1][2];
      
      u_xi = u_load[1] * (-1) + u_load[2] * (-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1] * (-1) + u[t_now][x][y][z][2] * (-1);
      
         
      f_e = 0.027778 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );
      
      f_e_o = 0.083333 * p_load * u_xi;
      
      f_x_t = 0.083333 * p[t_now][x][y][z] * u_x_xi;
      
      f[18] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;
}




