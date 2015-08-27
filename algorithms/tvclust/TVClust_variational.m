function [assignment] = TVClust_variational(C, SM, X)
dlmwrite('C2.txt',C,'delimiter',' ','precision',3)
dlmwrite('x2.txt',X,'delimiter',' ','precision',3)
dlmwrite('SM2.txt',SM,'delimiter',' ','precision',3)

system('Rscript main_VB_TVClust_sim_data2.r '); 
fileID = fopen('result2.txt','r');
assignment = fscanf(fileID,'%d'); 

end