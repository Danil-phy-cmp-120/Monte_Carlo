#include <iostream>
#include <fstream> //запись в файл
#include <ctime>
#include <cmath>
#include <cstdlib>
#include "omp.h"
#include <sstream>


using namespace std;
 
int main()
{
	
	float T_start, T_end, T_step;
    int nmcs, nstep, mcs_start, L, IQ_n, IQ_diff, J_size, J_size2;
    
    ofstream out_test;
    out_test.open("Test.dat");

    ifstream f1("input/Settings.dat");
    
    if (f1.is_open())//Если открытие файла прошло успешно
	{
		string line;
		
        while(getline(f1, line))
        {
            istringstream iss(line); //Создадим поток для считывания данных из строчки
            
            iss >> T_start >> T_end >> T_step >> nmcs >> nstep >> mcs_start >> L >> IQ_n >> IQ_diff >> J_size >> J_size2;
            //cout << T_start << "  " << T_end << "  " << T_step << "  " << nmcs << "  " << nstep << "  " << mcs_start << "  " << L << "  " << IQ_n << "  " << J_size << "\n";
            //out_test << T_start << "  " << T_end << "  " << T_step << "  " << nmcs << "  " << nstep << "  " << mcs_start << "  " << L << "  " << IQ_n << "  " << J_size << "\n";
            
        }
    }
    f1.close();
    
    //std::cout << J_size2 << "\n";
    
    float Lattice_cpp[L*L*L*IQ_n][7];
    float J_cpp[J_size][J_size2];
    int ind1, ind2;
    
    
    ind1 = 0;
    ind2 = 0;
    
    ifstream f2("input/J.dat");
   
	if (f2.is_open())
	{	
	    while(f2)
	    {
		    f2 >> J_cpp[ind1][ind2];
            
            ind2 = ind2 + 1;
            if (ind2 >= J_size2)
            {
				ind1 = ind1 + 1;
				ind2 = 0;
			}    
        }
    }
    f2.close();
        

    ind1 = 0;
    ind2 = 0;
    
    ifstream f3("input/Lattice.dat");
   
	if (f3.is_open())//Если открытие файла прошло успешно
	{	
	    while(f3)
	    {
		    f3 >> Lattice_cpp[ind1][ind2];
            
            ind2 = ind2 + 1;
            if (ind2> 6)
            {
				ind1 = ind1 + 1;
				ind2 = 0;
			}    
        }
    }
    f3.close();
   
    float magmom[IQ_diff];
    ind1 = 0;
    
    ifstream f4("input/magmom.dat");
	if (f4.is_open())
	{	
	    while(f4)
	    {
		    f4 >> magmom[ind1];
		    ind1 = ind1 + 1;   
        }
    }
    f4.close();
    
    int mri, mcs = 0, mcs_count = 0, num, conc[10];
	float E_sys = 0, kB, IQ_curent_atom;
	float E_1 = 0.0, E_2 = 0.0, E_calculate = 0.0;
	float Energy[int((T_end-T_start)/T_step) + 1][2];
	bool flag;
	float W, W1, W2, W_r, W_m, m, m1, m2, x0, y0, z0, x, y, z, H1 = 0.0, H2 = 0.0;
    float curentspin[3], spin_neighbors[J_size][3];
    float a[3], b[3], s[3];
    int num_neighbors = 0;
    float sum_spin_x = 0.0, sum_spin_y = 0.0, sum_spin_z = 0.0;
    float M_element[int((T_end-T_start)/T_step) + 1][IQ_diff*3], M_element_square[int((T_end-T_start)/T_step) + 1][IQ_diff*3], M_calculate[IQ_diff*3], M_all[int((T_end-T_start)/T_step) + 1][IQ_diff];
    float M_all_sum = 0.0;
    
    // Динамические массивы //  
    int **Neighbors_i = new int* [L*L*L*IQ_n]; // строки
    for (int count = 0; count < L*L*L*IQ_n; count++)
        Neighbors_i[count] = new int [J_size]; // столбцы
        
    int **Neighbors_j = new int* [L*L*L*IQ_n]; // строки
    for (int count = 0; count < L*L*L*IQ_n; count++)
        Neighbors_j[count] = new int [J_size]; // столбцы
    
    kB = 0.086173303; //мэВ / К
	

    // Обнуление массивов
    for ( int i = 0; i < J_size; i++) {
        for ( int j = 0; j < 3; j++) {
            spin_neighbors[i][j] = 0.0;
        }
    }
    for ( int i = 0; i < int((T_end-T_start)/T_step) + 1; i++) {
        for ( int j = 0; j < IQ_diff*3; j++) {
            M_element[i][j] = 0.0;
            M_element_square[i][j] = 0.0;
        }
    }
    for ( int i = 0; i < 3; i++) {
        curentspin[i] = 0.0;
    }
    
    for ( int i = 0; i < L*L*L*IQ_n; i++) {
        for ( int j = 0; j < J_size; j++) {
            Neighbors_i[i][j] = -1;
            Neighbors_j[i][j] = -1;
        }
    }
    

    #pragma omp parallel for firstprivate(Lattice_cpp, J_cpp) private(num_neighbors, IQ_curent_atom, x0, y0, z0, x, y, z, flag)
    for ( int k = 0; k < L*L*L*IQ_n; k++) {
		
		// Подсчет колическтва атомов каждого сорта
        conc[int(Lattice_cpp[k][0])] = conc[int(Lattice_cpp[k][0])] + 1;
    
	    //out_test << k << " / " << L*L*L*IQ_n << "\n";
	    std::cout << k << " / " << L*L*L*IQ_n << "\n";
        num_neighbors = 0;
    
        IQ_curent_atom = Lattice_cpp[k][0]; x0 = Lattice_cpp[k][1]; y0 = Lattice_cpp[k][2]; z0 = Lattice_cpp[k][3];
  
        for ( int i = 0; i < J_size; i++) {
  	        flag = false;
   	            		
   	        x = x0 + J_cpp[i][2];
            y = y0 + J_cpp[i][3];
   	        z = z0 + J_cpp[i][4];
   	        		    
   	        for ( int j = 0; j < L*L*L*IQ_n; j++) {
   	            if ((int(J_cpp[i][0]) == int(IQ_curent_atom) && int(J_cpp[i][1]) == int(Lattice_cpp[j][0]))) { //|| (int(J_cpp[i][1]) == int(IQ_curent_atom) && int(J_cpp[i][0]) == int(Lattice_cpp[j][0]))){
   	                if ((Lattice_cpp[j][1] == x) && (Lattice_cpp[j][2] == y) && (Lattice_cpp[j][3] == z)) {
					 	Neighbors_i[k][num_neighbors] = i;
					  	Neighbors_j[k][num_neighbors] = j;
                        num_neighbors = num_neighbors + 1;
               	        flag = true;
       	            }
   	            }  
   	        }
   	        		    
   	        if (flag == false){
	            for ( int q = 0; q < 25; q++) {
       	            for ( int j = 0; j < L*L*L*IQ_n; j++) {
						
						if (q == 0) {
    	  	                x = x0 + L + J_cpp[i][2];
                            y = y0 + J_cpp[i][3];
                            z = z0 + J_cpp[i][4];
					    }
				        if (q == 1) {
  	  	                    x = x0 + J_cpp[i][2];
                            y = y0 + L + J_cpp[i][3];
                            z = z0 + J_cpp[i][4];
				        }
    			        if (q == 2) {
   	  	                    x = x0 + J_cpp[i][2];
                            y = y0 + J_cpp[i][3];
                            z = z0 + L + J_cpp[i][4];
				        }
				        if (q == 3) {
  	  	                    x = x0 + L + J_cpp[i][2];
                            y = y0 + L + J_cpp[i][3];
                            z = z0 + J_cpp[i][4];
				        }
				        if (q == 4) {
   	  	                    x = x0 + L + J_cpp[i][2];
                            y = y0 + J_cpp[i][3];
                            z = z0 + L + J_cpp[i][4];
				        }
				        if (q == 5) {
   	  	                    x = x0 + J_cpp[i][2];
                            y = y0 + L + J_cpp[i][3];
                            z = z0 + L + J_cpp[i][4];
				        }
				        if (q == 6) {
  	  	                    x = x0 + L + J_cpp[i][2];
                            y = y0 + L + J_cpp[i][3];
                            z = z0 + L + J_cpp[i][4];
				        }
    			        if (q == 7) {
   	  	                    x = x0 - L + J_cpp[i][2];
                            y = y0 + J_cpp[i][3];
                            z = z0 + J_cpp[i][4];
				        }
				        if (q == 8) {
  	  	                    x = x0 + J_cpp[i][2];
                            y = y0 - L + J_cpp[i][3];
                            z = z0 + J_cpp[i][4];
				        }
				        if (q == 9) {
   	  	                    x = x0 + J_cpp[i][2];
                            y = y0 + J_cpp[i][3];
                            z = z0 - L + J_cpp[i][4];
				        }
				        if (q == 10) {
  	  	                    x = x0 - L + J_cpp[i][2];
                            y = y0 - L + J_cpp[i][3];
                            z = z0 + J_cpp[i][4];
				        }
				        if (q == 11) {
  	  	                    x = x0 - L + J_cpp[i][2];
                            y = y0 + J_cpp[i][3];
                            z = z0 - L + J_cpp[i][4];
				        }
				        if (q == 12) {
   	  	                    x = x0 + J_cpp[i][2];
                            y = y0 - L + J_cpp[i][3];
                            z = z0 - L + J_cpp[i][4];
				        }
				        if (q == 13) {
  	  	                    x = x0 - L + J_cpp[i][2];
                            y = y0 - L + J_cpp[i][3];
                            z = z0 - L + J_cpp[i][4];
				        }
				        if (q == 14) {
  	  	                    x = x0 + L + J_cpp[i][2];
                            y = y0 - L + J_cpp[i][3];
                            z = z0 + J_cpp[i][4];
                        }
                        if (q == 15) {
     	  	                x = x0 - L + J_cpp[i][2];
                            y = y0 + L + J_cpp[i][3];
                            z = z0 + J_cpp[i][4];
                        }
                        if (q == 16) {
     	                    x = x0 - L + J_cpp[i][2];
                            y = y0 + J_cpp[i][3];
                            z = z0 + L + J_cpp[i][4];
						}
                        if (q == 17) {
    	                    x = x0 + L + J_cpp[i][2];
                            y = y0 + J_cpp[i][3];
                            z = z0 - L + J_cpp[i][4];
						}
                        if (q == 18) {
    	                    x = x0 + J_cpp[i][2];
                            y = y0 + L + J_cpp[i][3];
                            z = z0 - L + J_cpp[i][4];
						}
                        if (q == 19) {
      	                    x = x0 + J_cpp[i][2];
                            y = y0 - L + J_cpp[i][3];
                            z = z0 + L + J_cpp[i][4];
						}
                        if (q == 20) {
 	  	                    x = x0 - L + J_cpp[i][2];
                            y = y0 - L + J_cpp[i][3];
                            z = z0 + L + J_cpp[i][4];
						}
                        if (q == 21) {
 	  	                    x = x0 - L + J_cpp[i][2];
                            y = y0 + L + J_cpp[i][3];
                            z = z0 - L + J_cpp[i][4];
						}
                        if (q == 22) {
 	  	                    x = x0 + L + J_cpp[i][2];
                            y = y0 - L + J_cpp[i][3];
                            z = z0 - L + J_cpp[i][4];
						}
                        if (q == 23) {
   	  	                    x = x0 + L + J_cpp[i][2];
                            y = y0 + L + J_cpp[i][3];
                            z = z0 - L + J_cpp[i][4];
						}
                        if (q == 24) {
  	  	                    x = x0 + L + J_cpp[i][2];
                            y = y0 - L + J_cpp[i][3];
                            z = z0 + L + J_cpp[i][4];
						}
                        if (q == 25) {
   	  	                    x = x0 - L + J_cpp[i][2];
                            y = y0 + L + J_cpp[i][3];
                            z = z0 + L + J_cpp[i][4];
						}

			            if ((int(J_cpp[i][0]) == int(IQ_curent_atom) && int(J_cpp[i][1]) == int(Lattice_cpp[j][0]))) {  // || (int(J_cpp[i][1]) == int(IQ_curent_atom) && int(J_cpp[i][0]) == int(Lattice_cpp[j][0]))){
  		                    if ((Lattice_cpp[j][1] == x) && (Lattice_cpp[j][2] == y) && (Lattice_cpp[j][3] == z)) {
                                Neighbors_i[k][num_neighbors] = i;
     				            Neighbors_j[k][num_neighbors] = j;
					            //out << "Welcome to CPP" << std::endl;
                                num_neighbors = num_neighbors + 1;
                            }
                        }		
	                }
		        }
		    }
		}
	
	    ofstream out1;
	    out1.open("input/Neighbors_i.dat");
        ofstream out2;
	    out2.open("input/Neighbors_j.dat");
	    for ( int i = 0; i < L*L*L*IQ_n; i++) {
            for ( int j = 0; j < J_size; j++) {
				if (Neighbors_j[i][j] != -1) {
                    out1 << Neighbors_i[i][j] << "  ";
                    out2 << Neighbors_j[i][j] << "  ";
                }
            }
            out1 << "\n";
            out2 << "\n";
        }
        out1.close(); 
        out2.close();
    }
	    
	srand((unsigned)time(NULL)); // автоматическая рандомизация     
    #pragma omp parallel for firstprivate(Lattice_cpp, J_cpp, curentspin, spin_neighbors, Neighbors_i, Neighbors_j, M_calculate) private(flag, mcs_count, mcs, W, mri, IQ_curent_atom, x0, y0, z0, x, y, z, H1, H2, W1, W2, m1, m2, m, W_r, W_m, E_calculate, E_1, E_2, num_neighbors)
    
    for (int t = int(T_start); t <= int(T_end); t = t + int(T_step)) {
		mcs_count = 0; mcs = 0;
		for ( int n = 1; n <= nmcs; n++) { 		
			//out_test << t << "   " << n << "   " << omp_get_thread_num() << "\n";
			out_test << t << "   " << n << "\n";
			//std::cout << t << "   " << n << "   " << omp_get_thread_num() << "\n";
			//std::cout << n << "   "  << "\n";
    		for ( int count = 1; count <= L*L*L*IQ_n; count++) { //Один шаг Монте-Карло
    			
    			W = (double)rand() / RAND_MAX; mri = int(L*L*L*(IQ_n)*W); 
    			if (mri == L*L*L*IQ_n) {
					mri = mri - 1;
				}
	          
	    		for ( int i = 0; i < 3; i++) {
		    	    curentspin[i] = Lattice_cpp[mri][i+4];
				}	
    			
    			for ( int num_j = 0; num_j < J_size; num_j++) {
					if (Neighbors_j[mri][num_j] != -1) {
						spin_neighbors[Neighbors_i[mri][num_j]][0] = Lattice_cpp[Neighbors_j[mri][num_j]][4]; 
           	            spin_neighbors[Neighbors_i[mri][num_j]][1] = Lattice_cpp[Neighbors_j[mri][num_j]][5];
                        spin_neighbors[Neighbors_i[mri][num_j]][2] = Lattice_cpp[Neighbors_j[mri][num_j]][6];
					} 
				}

	            for ( int num_j = 0; num_j < J_size; num_j++) {
				    if (Neighbors_j[mri][num_j] != -1) {
						if (J_size2 <= 7) {
							H1 = H1 - (J_cpp[Neighbors_i[mri][num_j]][6]/kB) * (spin_neighbors[Neighbors_i[mri][num_j]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][num_j]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][num_j]][2] * curentspin[2]);
					    }
					    else {
							H1 = H1 - (J_cpp[Neighbors_i[mri][num_j]][6]*t*t + J_cpp[Neighbors_i[mri][num_j]][7]*t + J_cpp[Neighbors_i[mri][num_j]][8])/kB * (spin_neighbors[Neighbors_i[mri][num_j]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][num_j]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][num_j]][2] * curentspin[2]);		
						}
					    //out_test << Lattice_cpp[mri][0] << "  " << Lattice_cpp[Neighbors_j[mri][num_j]][0] << "   " << (J_cpp[Neighbors_i[mri][num_j]][6]/kB) << "   " << spin_neighbors[Neighbors_i[mri][num_j]][0] << "   " << spin_neighbors[Neighbors_i[mri][num_j]][1] << "   " << spin_neighbors[Neighbors_i[mri][num_j]][2] << "   " << curentspin[0] << "   " << curentspin[1] << "   " << curentspin[2] << "\n";
					}	
   				}
   				//out_test << H1 << "\n" << "\n";
   				
   				//Изинг
        		//curentspin[2] = -curentspin[2];
        
    			//Гейзенберг            
				W1 = (double)rand() / RAND_MAX; 
				W2 = (double)rand() / RAND_MAX; 
				m1 = 1.0 - 2.0 * W1;
    			m2 = 1.0 - 2.0 * W2;
    			m = sqrt(m1*m1 + m2*m2); 
    	
    			if (m*m < 1){
        	        curentspin[0] = 2.0 * m1 * sqrt(1.0 - m*m);
        	        curentspin[1] = 2.0 * m2 * sqrt(1.0 - m*m);
        		    curentspin[2] = 1.0 - 2.0*m*m;
    			} 
    	       
	            for ( int num_j = 0; num_j < J_size; num_j++) {
				    if (Neighbors_j[mri][num_j] != -1) {
					   	if (J_size2 <= 7) {
							H2 = H2 - (J_cpp[Neighbors_i[mri][num_j]][6]/kB) * (spin_neighbors[Neighbors_i[mri][num_j]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][num_j]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][num_j]][2] * curentspin[2]);		
					    }
					    else {
							H2 = H2 - (J_cpp[Neighbors_i[mri][num_j]][6]*t*t + J_cpp[Neighbors_i[mri][num_j]][7]*t + J_cpp[Neighbors_i[mri][num_j]][8])/kB * (spin_neighbors[Neighbors_i[mri][num_j]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][num_j]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][num_j]][2] * curentspin[2]);		
						}
					    //out_test << Lattice_cpp[mri][0] << "  " << Lattice_cpp[Neighbors_j[mri][num_j]][0] << "   " << (J_cpp[Neighbors_i[mri][num_j]][6]/kB) << "   " << spin_neighbors[Neighbors_i[mri][num_j]][0] << "   " << spin_neighbors[Neighbors_i[mri][num_j]][1] << "   " << spin_neighbors[Neighbors_i[mri][num_j]][2] << "   " << curentspin[0] << "   " << curentspin[1] << "   " << curentspin[2] << "\n";
    			    }
			    }
    		
    	        //out_test << count << "  " << Lattice_cpp[mri][0] << "  " << H1 << "   " << H2  << "\n";  
    	        //std::cout << H1 << "   " << H2 << '\n';
    	        
    	        
    			if ((H2 - H1) < 0){  
					Lattice_cpp[mri][4] = curentspin[0];
					Lattice_cpp[mri][5] = curentspin[1];
					Lattice_cpp[mri][6] = curentspin[2];
					//out_test << mri << Lattice_cpp[mri][4] << "   " << Lattice_cpp[mri][5] << "   " << Lattice_cpp[mri][6] << "\n"<< "\n";
				}
				else {
	    			W_r = (double)rand() / RAND_MAX;	
	    			W_m = exp(-2*(H2 - H1)/t);
	    			if (W_r < W_m){  
	        			Lattice_cpp[mri][4] = curentspin[0];
		    			Lattice_cpp[mri][5] = curentspin[1];
			    		Lattice_cpp[mri][6] = curentspin[2];	
		    		}
		    		//out_test << "else" << Lattice_cpp[mri][4] << "   " << Lattice_cpp[mri][5] << "   " << Lattice_cpp[mri][6] << "\n"<< "\n";
				}  
				H1 = 0.0;
				H2 = 0.0;
				
				for ( int i = 0; i < J_size; i++) {
                    for ( int j = 0; j < 3; j++) {
                        spin_neighbors[i][j] = 0.0;
                    }
                }
                
			}
			//out_test << "\n";
			

    		mcs = mcs + 1; 
    		if (mcs >= mcs_start) mcs_count = mcs_count + 1;
             
            if (mcs_count == nstep){ 
			    mcs_count = 0;	
			    E_calculate	 = 0.0;
			    for ( int i = 0; i < IQ_diff*3; i++) {
		    	    M_calculate[i] = 0.0;
				}
			    
                ////////////////////////////////////////////////// Рассчет E /////////////////////////////////////////////////////////////////
		        for ( int k = 0; k < L*L*L*IQ_n; k++) {
    
	    		    for ( int i = 0; i < 3; i++) {
		    	        curentspin[i] = Lattice_cpp[k][i+4];
				    }
         
    			    for ( int num_j = 0; num_j < J_size; num_j++) {
					    if (Neighbors_j[k][num_j] != -1) {
						    spin_neighbors[Neighbors_i[k][num_j]][0] = Lattice_cpp[Neighbors_j[k][num_j]][4]; 
           	                spin_neighbors[Neighbors_i[k][num_j]][1] = Lattice_cpp[Neighbors_j[k][num_j]][5];
                            spin_neighbors[Neighbors_i[k][num_j]][2] = Lattice_cpp[Neighbors_j[k][num_j]][6];
					    }
				    }
				    
				    H1 = 0;
	                for ( int num_j = 0; num_j < J_size; num_j++) {
				        if (Neighbors_j[k][num_j] != -1) {
					        if (J_size2 <= 7) {
							    H1 = H1 - (J_cpp[Neighbors_i[mri][num_j]][6]/kB) * (spin_neighbors[Neighbors_i[mri][num_j]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][num_j]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][num_j]][2] * curentspin[2]);		
					        }
					        else {
						    	H1 = H1 - (J_cpp[Neighbors_i[mri][num_j]][6]*t*t + J_cpp[Neighbors_i[mri][num_j]][7]*t + J_cpp[Neighbors_i[mri][num_j]][8])/kB * (spin_neighbors[Neighbors_i[mri][num_j]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][num_j]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][num_j]][2] * curentspin[2]);		
						    }
					        //out_test << Lattice_cpp[k][0] << "  " << Lattice_cpp[Neighbors_j[k][num_j]][0] << "   " << (J_cpp[Neighbors_i[k][num_j]][6]/kB) << "   " << spin_neighbors[Neighbors_i[k][num_j]][0] << "   " << spin_neighbors[Neighbors_i[k][num_j]][1] << "   " << spin_neighbors[Neighbors_i[k][num_j]][2] << "   " << curentspin[0] << "   " << curentspin[1] << "   " << curentspin[2] << "\n";
					    }	
   				    }
   			        E_calculate = E_calculate + 0.5*H1; 	 
   			        
   			        ////////////////////////////////// Рассчет M ////////////////////////////////////////////////////
		    	    for ( int i = 0; i < IQ_diff; i++) {
		    	        if (Lattice_cpp[k][0] == (i+1)) {
		    	            M_calculate[3*i] = M_calculate[3*i] + Lattice_cpp[k][4];
		    	            M_calculate[1+3*i] = M_calculate[1+3*i] + Lattice_cpp[k][5];
		    	            M_calculate[2+3*i] = M_calculate[2+3*i] + Lattice_cpp[k][6];
					    }
				    }
	            }
	            
		        E_1 = E_1 + E_calculate;					
                E_2 = E_2 + E_calculate*E_calculate;

                for ( int i = 0; i < IQ_diff*3; i++) {
					M_element[int((t-T_start)/T_step)][i] = M_element[int((t-T_start)/T_step)][i] + M_calculate[i];
					M_element_square[int((t-T_start)/T_step)][i] = M_element_square[int((t-T_start)/T_step)][i] + M_calculate[i]*M_calculate[i];
				}
                
    	    }
	    }
    	
        E_1 = E_1/(L*L*L*IQ_n);
        E_2 = E_2/(L*L*L*IQ_n);
	                    
	    E_1 = E_1 /((nmcs - mcs_start)/nstep);
	    E_2 = E_2 /((nmcs - mcs_start)/nstep);    
	    
        for ( int j = 0; j < IQ_diff*3; j++) {
			    M_element[int((t-T_start)/T_step)][j] = M_element[int((t-T_start)/T_step)][j]/(L*L*L);
			    M_element_square[int((t-T_start)/T_step)][j] = M_element_square[int((t-T_start)/T_step)][j]/(L*L*L);
			    M_element[int((t-T_start)/T_step)][j] = M_element[int((t-T_start)/T_step)][j]/((nmcs - mcs_start)/nstep);
			    M_element_square[int((t-T_start)/T_step)][j] = M_element_square[int((t-T_start)/T_step)][j]/((nmcs - mcs_start)/nstep);
        }
	    
        Energy[int((t-T_start)/T_step)][0] = E_1; 
        Energy[int((t-T_start)/T_step)][1] = E_2; 
    
        E_1 = 0.0; E_2 = 0.0;
        
        // Расчет полной намагниченности //
        for ( int i = 0; i < IQ_diff; i++) {
			for ( int j = 0; j < 3; j++) {
                M_all[int((t-T_start)/T_step)][i] = M_all[int((t-T_start)/T_step)][i] + M_element[int((t-T_start)/T_step)][j+3*i]*M_element[int((t-T_start)/T_step)][j+3*i];
			}
			M_all[int((t-T_start)/T_step)][i] = sqrt(M_all[int((t-T_start)/T_step)][i]);
		}
    }

    // Вывод результатов в файлы //
    ofstream out_E;
	out_E.open("Energy.dat");
	for ( int i = 0; i < int((T_end-T_start)/T_step) + 1; i++) {
		out_E << T_start + i*T_step  << "\t";
	    for ( int j = 0; j < 2; j++) {
            out_E << Energy[i][j] << "\t";
        }
        out_E << "\n";
    }
    out_E.close();
    
    ofstream out_M;
	out_M.open("M_element.dat");
	for ( int i = 0; i < int((T_end-T_start)/T_step) + 1; i++) {
		out_M << T_start + i*T_step  << "\t";
	    for ( int j = 0; j < IQ_diff*3; j++) {
            out_M << M_element[i][j] << "\t";
        }
        out_M << "\n";
    }
    out_M.close();

    ofstream out_M2;
	out_M2.open("M_element_square.dat");
	for ( int i = 0; i < int((T_end-T_start)/T_step) + 1; i++) {
		out_M2 << T_start + i*T_step  << "\t";
	    for ( int j = 0; j < IQ_diff*3; j++) {
            out_M2 << M_element_square[i][j] << "\t";
        }
        out_M2 << "\n";
    }
    out_M2.close();
    
    ofstream out_M_all;
	out_M_all.open("M_all.dat");
	for ( int i = 0; i < int((T_end-T_start)/T_step) + 1; i++) {
		M_all_sum = 0.0;
		out_M_all << T_start + i*T_step << "\t";
		for ( int j = 0; j < IQ_diff; j++) {
			M_all_sum = M_all_sum + M_all[i][j] * magmom[j];
		}
		out_M_all << M_all_sum  << "\n";
    }
    out_M_all.close();
    
    out_test.close();
    
    for (int count = 0; count < L*L*L*IQ_n; count++) {
        delete [] Neighbors_i[count];
        delete [] Neighbors_j[count];
	}
	delete [] Neighbors_i;
    delete [] Neighbors_j;
    
}
