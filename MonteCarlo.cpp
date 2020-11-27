#include <iostream>
#include <fstream> //запись в файл
#include <ctime>
#include <cmath>
#include <cstdlib>
#include "omp.h"
#include <sstream>
#include <random>
#include <vector>


using namespace std;
 
int main()
{
	ofstream out_test;
    out_test.open("Test.dat");
	
	
	float T_start, T_end, T_step;
    int nmcs, nstep, mcs_start, L, IQ_n, IQ_diff, J_size, J_size2;
    ifstream f1("input/Settings.dat");
    if (f1.is_open())
	{
		string line;
        while(getline(f1, line))
        {
            istringstream iss(line);
            iss >> T_start >> T_end >> T_step >> nmcs >> nstep >> mcs_start >> L >> IQ_n >> IQ_diff >> J_size >> J_size2;       
        }
    }
    f1.close();
    
    
    std::vector<std::vector<float> > J_cpp(J_size, std::vector<float>(J_size2));
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
        

    std::vector<std::vector<float> > Lattice_cpp(L*L*L*IQ_n, std::vector<float>(7));    
    
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

    std::vector<float> magmom = std::vector<float>(IQ_diff);
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
    
    
    
    
    
    int mri, mcs = 0, mcs_count = 0, num, iterator, num_neighbors = 0;
	float E_sys = 0, kB = 0.086173303, IQ_curent_atom, E_1 = 0.0, E_2 = 0.0, E_calculate = 0.0, W, W1, W2, W_r, W_m, m, m1, m2, x0, y0, z0, x, y, z, H1 = 0.0, H2 = 0.0, M_all_sum = 0.0, curentspin[3];
	bool flag;
	
	// Вектора //
	std::vector<std::vector<int> > Neighbors_i(L*L*L*IQ_n, std::vector<int>(J_size));
	std::vector<std::vector<int> > Neighbors_j(L*L*L*IQ_n, std::vector<int>(J_size));
	std::vector<std::vector<float> > M_all(int((T_end-T_start)/T_step) + 1, std::vector<float>(IQ_diff));
	std::vector<std::vector<float> > M_element(int((T_end-T_start)/T_step) + 1, std::vector<float>(IQ_diff*3));
	std::vector<std::vector<float> > M_element_square(int((T_end-T_start)/T_step) + 1, std::vector<float>(IQ_diff*3));
	std::vector<std::vector<float> > Energy(int((T_end-T_start)/T_step) + 1, std::vector<float>(2));
    std::vector<std::vector<float> > spin_neighbors(J_size, std::vector<float>(3));
    std::vector<float> M_calculate = std::vector<float>(IQ_diff*3);
 
    /*// Динамические массивы //  
    int **Neighbors_i = new int* [L*L*L*IQ_n]; 
    for (int count = 0; count < L*L*L*IQ_n; count++)
        Neighbors_i[count] = new int [J_size]; 
        
    int **Neighbors_j = new int* [L*L*L*IQ_n]; 
    for (int count = 0; count < L*L*L*IQ_n; count++)
        Neighbors_j[count] = new int [J_size]; 
        
    float **M_all = new float* [int((T_end-T_start)/T_step) + 1]; 
    for (int count = 0; count < int((T_end-T_start)/T_step) + 1; count++)
        M_all[count] = new float [IQ_diff];
    
    float **M_element = new float* [int((T_end-T_start)/T_step) + 1]; 
    for (int count = 0; count < int((T_end-T_start)/T_step) + 1; count++)
        M_element[count] = new float [IQ_diff*3]; 
        
    float **M_element_square = new float* [int((T_end-T_start)/T_step) + 1]; 
    for (int count = 0; count < int((T_end-T_start)/T_step) + 1; count++)
        M_element_square[count] = new float [IQ_diff*3]; 
        
    float **Energy = new float* [int((T_end-T_start)/T_step) + 1];
    for (int count = 0; count < int((T_end-T_start)/T_step) + 1; count++)
        Energy[count] = new float [2]; 
        
    float **spin_neighbors = new float* [J_size];
    for (int count = 0; count < J_size; count++)
        spin_neighbors[count] = new float [3]; 
        
    float *M_calculate = new float [IQ_diff*3];*/


    // Обнуление массивов

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
    
    for ( int i = 0; i < J_size; i++) {
        for ( int j = 0; j < 3; j++) {
            spin_neighbors[i][j] = 0.0;
        }
    }
    

    #pragma omp parallel for private(num_neighbors, IQ_curent_atom, x0, y0, z0, x, y, z, flag)
    for ( int k = 0; k < L*L*L*IQ_n; k++) {
    
	    //out_test << k << " / " << L*L*L*IQ_n << "\n";
	    std::cout << k << " / " << L*L*L*IQ_n << "\n";
        num_neighbors = 0;
    
        IQ_curent_atom = Lattice_cpp[k][0]; x0 = Lattice_cpp[k][1]; y0 = Lattice_cpp[k][2]; z0 = Lattice_cpp[k][3];
  
        for ( int i = 0; i < J_size; i++) {
  	        flag = false;
   	            		
   	        x = x0 + J_cpp[i][2]; y = y0 + J_cpp[i][3]; z = z0 + J_cpp[i][4];
   	        		    
   	        for ( int j = 0; j < L*L*L*IQ_n; j++) {
   	            if ((int(J_cpp[i][0]) == int(IQ_curent_atom) && int(J_cpp[i][1]) == int(Lattice_cpp[j][0]))) {
   	                if ((Lattice_cpp[j][1] == x) && (Lattice_cpp[j][2] == y) && (Lattice_cpp[j][3] == z)) {
					 	Neighbors_i[k][num_neighbors] = i;
					  	Neighbors_j[k][num_neighbors] = j;
                        num_neighbors = num_neighbors + 1;
               	        flag = true;
       	            }
   	            }  
   	        }
   	        		    
   	        if (flag == false){
	            for ( int q = 0; q <= 25; q++) {
       	            for ( int j = 0; j < L*L*L*IQ_n; j++) {
				        switch (q)
                        {
						    case 0: x = x0 + L + J_cpp[i][2]; y = y0 + J_cpp[i][3]; z = z0 + J_cpp[i][4]; break;
				            case 1: x = x0 + J_cpp[i][2]; y = y0 + L + J_cpp[i][3]; z = z0 + J_cpp[i][4]; break;
				            case 2: x = x0 + J_cpp[i][2]; y = y0 + J_cpp[i][3]; z = z0 + L + J_cpp[i][4]; break;
                            case 3: x = x0 + L + J_cpp[i][2]; y = y0 + L + J_cpp[i][3]; z = z0 + J_cpp[i][4]; break;
				            case 4: x = x0 + L + J_cpp[i][2]; y = y0 + J_cpp[i][3]; z = z0 + L + J_cpp[i][4]; break;
				            case 5: x = x0 + J_cpp[i][2]; y = y0 + L + J_cpp[i][3]; z = z0 + L + J_cpp[i][4]; break;
				            case 6: x = x0 + L + J_cpp[i][2]; y = y0 + L + J_cpp[i][3]; z = z0 + L + J_cpp[i][4]; break;
				            case 7: x = x0 - L + J_cpp[i][2]; y = y0 + J_cpp[i][3]; z = z0 + J_cpp[i][4]; break;
				            case 8: x = x0 + J_cpp[i][2]; y = y0 - L + J_cpp[i][3]; z = z0 + J_cpp[i][4]; break;
				            case 9: x = x0 + J_cpp[i][2]; y = y0 + J_cpp[i][3]; z = z0 - L + J_cpp[i][4]; break;
				            case 10: x = x0 - L + J_cpp[i][2]; y = y0 - L + J_cpp[i][3]; z = z0 + J_cpp[i][4]; break;
				            case 11: x = x0 - L + J_cpp[i][2]; y = y0 + J_cpp[i][3]; z = z0 - L + J_cpp[i][4]; break;
				            case 12: x = x0 + J_cpp[i][2]; y = y0 - L + J_cpp[i][3]; z = z0 - L + J_cpp[i][4]; break;
				            case 13: x = x0 - L + J_cpp[i][2]; y = y0 - L + J_cpp[i][3]; z = z0 - L + J_cpp[i][4]; break;
				            case 14: x = x0 + L + J_cpp[i][2]; y = y0 - L + J_cpp[i][3]; z = z0 + J_cpp[i][4]; break;
                            case 15: x = x0 - L + J_cpp[i][2]; y = y0 + L + J_cpp[i][3]; z = z0 + J_cpp[i][4]; break;
                            case 16: x = x0 - L + J_cpp[i][2]; y = y0 + J_cpp[i][3]; z = z0 + L + J_cpp[i][4]; break;
					        case 17: x = x0 + L + J_cpp[i][2]; y = y0 + J_cpp[i][3]; z = z0 - L + J_cpp[i][4]; break;
					        case 18: x = x0 + J_cpp[i][2]; y = y0 + L + J_cpp[i][3]; z = z0 - L + J_cpp[i][4]; break;
					        case 19: x = x0 + J_cpp[i][2]; y = y0 - L + J_cpp[i][3]; z = z0 + L + J_cpp[i][4]; break;
					        case 20: x = x0 - L + J_cpp[i][2]; y = y0 - L + J_cpp[i][3]; z = z0 + L + J_cpp[i][4]; break;
					        case 21: x = x0 - L + J_cpp[i][2]; y = y0 + L + J_cpp[i][3]; z = z0 - L + J_cpp[i][4]; break;
					        case 22: x = x0 + L + J_cpp[i][2]; y = y0 - L + J_cpp[i][3]; z = z0 - L + J_cpp[i][4]; break;
					        case 23: x = x0 + L + J_cpp[i][2]; y = y0 + L + J_cpp[i][3]; z = z0 - L + J_cpp[i][4]; break;
					        case 24: x = x0 + L + J_cpp[i][2]; y = y0 - L + J_cpp[i][3]; z = z0 + L + J_cpp[i][4]; break;
					        case 25: x = x0 - L + J_cpp[i][2]; y = y0 + L + J_cpp[i][3]; z = z0 + L + J_cpp[i][4]; break;
					    }

			            if ((int(J_cpp[i][0]) == int(IQ_curent_atom) && int(J_cpp[i][1]) == int(Lattice_cpp[j][0]))) {
  		                    if ((Lattice_cpp[j][1] == x) && (Lattice_cpp[j][2] == y) && (Lattice_cpp[j][3] == z)) {
                                Neighbors_i[k][num_neighbors] = i;
     				            Neighbors_j[k][num_neighbors] = j;
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
	    
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, L*L*L*IQ_n-1); //Для выбора атома
    std::uniform_real_distribution<> dis2(0, 1); //Для выбора спина


    #pragma omp parallel for firstprivate(curentspin, Lattice_cpp, spin_neighbors) private(mcs_count, mcs, W, mri, H1, H2, W1, W2, m1, m2, m, W_r, W_m, E_calculate, E_1, E_2)
    for (int t = int(T_start); t <= int(T_end); t = t + int(T_step)) {
		mcs_count = 0; mcs = 0;
		for ( int n = 1; n <= nmcs; n++) { 		
			//out_test << t << "   " << n << "   " << omp_get_thread_num() << "\n";
			//out_test << t << "   " << n << "\n";
			std::cout << t << "   " << n << "   " << omp_get_thread_num() << "\n";
		    //std::cout << t << "   " << n << "\n";
    		for ( int count = 1; count <= L*L*L*IQ_n; count++) { //Один шаг Монте-Карло
				mri = dis(gen);
	            
	    		for ( int i = 0; i < 3; i++) {
		    	    curentspin[i] = Lattice_cpp[mri][i+4];
				}	
    			
    			for ( int i = 0; i < J_size; i++) {
			        if (Neighbors_j[mri][i] != -1) {
				        spin_neighbors[Neighbors_i[mri][i]][0] = Lattice_cpp[Neighbors_j[mri][i]][4]; 
           	            spin_neighbors[Neighbors_i[mri][i]][1] = Lattice_cpp[Neighbors_j[mri][i]][5];
                        spin_neighbors[Neighbors_i[mri][i]][2] = Lattice_cpp[Neighbors_j[mri][i]][6];
                    }
				}

	            for ( int i = 0; i < J_size; i++) {
				    if (Neighbors_j[mri][i] != -1) {
				        if (J_size2 <= 7) {
					        H1 = H1 - (J_cpp[Neighbors_i[mri][i]][6]/kB) * (spin_neighbors[Neighbors_i[mri][i]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][i]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][i]][2] * curentspin[2]);
					    }
					    else {
				            H1 = H1 - (J_cpp[Neighbors_i[mri][i]][6]*t*t + J_cpp[Neighbors_i[mri][i]][7]*t + J_cpp[Neighbors_i[mri][i]][8])/kB * (spin_neighbors[Neighbors_i[mri][i]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][i]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][i]][2] * curentspin[2]);		
				        }
				    }
                    else {
                        continue;
					}
   				}
   				
   				//Изинг
        		//curentspin[2] = -curentspin[2];
        
    			//Гейзенберг            
				W1 = dis2(gen);
				W2 = dis2(gen); 
				m1 = 1.0 - 2.0 * W1;
    			m2 = 1.0 - 2.0 * W2;
    			m = sqrt(m1*m1 + m2*m2); 
    	
    			if (m*m < 1){
        	        curentspin[0] = 2.0 * m1 * sqrt(1.0 - m*m);
        	        curentspin[1] = 2.0 * m2 * sqrt(1.0 - m*m);
        		    curentspin[2] = 1.0 - 2.0*m*m;
    			} 
    	       
	            for ( int i = 0; i < J_size; i++) {
				    if (Neighbors_j[mri][i] != -1) {
				        if (J_size2 <= 7) {
					        H2 = H2 - (J_cpp[Neighbors_i[mri][i]][6]/kB) * (spin_neighbors[Neighbors_i[mri][i]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][i]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][i]][2] * curentspin[2]);
					    }
					    else {
				            H2 = H2 - (J_cpp[Neighbors_i[mri][i]][6]*t*t + J_cpp[Neighbors_i[mri][i]][7]*t + J_cpp[Neighbors_i[mri][i]][8])/kB * (spin_neighbors[Neighbors_i[mri][i]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][i]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][i]][2] * curentspin[2]);		
				        }
				    }
                    else {
                        continue;
					}
   				}
    		
    	                
    			if ((H2 - H1) < 0){  
					Lattice_cpp[mri][4] = curentspin[0];
					Lattice_cpp[mri][5] = curentspin[1];
					Lattice_cpp[mri][6] = curentspin[2];
				}
				else {
	    			W_r = dis2(gen);
	    			W_m = exp(-2*(H2 - H1)/t);
	    			if (W_r < W_m){  
	        			Lattice_cpp[mri][4] = curentspin[0];
		    			Lattice_cpp[mri][5] = curentspin[1];
			    		Lattice_cpp[mri][6] = curentspin[2];	
		    		}
				}  
				H1 = 0.0;
				H2 = 0.0;
				
				for ( int i = 0; i < J_size; i++) {
                    for ( int j = 0; j < 3; j++) {
                        spin_neighbors[i][j] = 0.0;
                    }
                }
                
			}
			

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
         
    			    for ( int i = 0; i < J_size; i++) {
					    if (Neighbors_j[k][i] != -1) {
						    spin_neighbors[Neighbors_i[k][i]][0] = Lattice_cpp[Neighbors_j[k][i]][4]; 
           	                spin_neighbors[Neighbors_i[k][i]][1] = Lattice_cpp[Neighbors_j[k][i]][5];
                            spin_neighbors[Neighbors_i[k][i]][2] = Lattice_cpp[Neighbors_j[k][i]][6];
					    }
				    }
				    
				    H1 = 0;
	                for ( int i = 0; i < J_size; i++) {
				        if (Neighbors_j[mri][i] != -1) {
				            if (J_size2 <= 7) {
					            H1 = H1 - (J_cpp[Neighbors_i[mri][i]][6]/kB) * (spin_neighbors[Neighbors_i[mri][i]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][i]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][i]][2] * curentspin[2]);
					         }
					        else {
				                H1 = H1 - (J_cpp[Neighbors_i[mri][i]][6]*t*t + J_cpp[Neighbors_i[mri][i]][7]*t + J_cpp[Neighbors_i[mri][i]][8])/kB * (spin_neighbors[Neighbors_i[mri][i]][0] * curentspin[0] + spin_neighbors[Neighbors_i[mri][i]][1] * curentspin[1] + spin_neighbors[Neighbors_i[mri][i]][2] * curentspin[2]);		
				            }
						}
                        else {
                             continue;
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

    
    /*for (int count = 0; count < L*L*L*IQ_n; count++) {
        delete [] Neighbors_i[count];
        delete [] Neighbors_j[count];
	}
	delete [] Neighbors_i;
    delete [] Neighbors_j;
    
    for (int count = 0; count < int((T_end-T_start)/T_step) + 1; count++) {
	    delete [] Energy[count];
	    delete [] M_element[count];
	    delete [] M_element_square[count];
	    delete [] M_all[count];
    }
	delete [] Energy;
    delete [] M_element;
	delete [] M_element_square;
    delete [] M_all; 
    delete [] M_calculate; 
    
    for (int count = 0; count < J_size; count++) {
	    delete [] spin_neighbors[count];
    }
    delete [] spin_neighbors;     
    
    for (int count = 0; count < L*L*L*IQ_n; count++) {
		delete [] Lattice_cpp[count];
	    //delete [] Lattice_cpp_private[count];
	}
    delete [] Lattice_cpp;     
    //delete [] Lattice_cpp_private;  
    
    for (int count = 0; count < J_size; count++) {
		delete [] J_cpp[count];
	}
    delete [] J_cpp;  */      
}
