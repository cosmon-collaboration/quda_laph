with(LinearAlgebra):

sigma1:=Matrix([[0,1],[1,0]]);
sigma2:=Matrix([[0,-I],[I,0]]);
sigma3:=Matrix([[1,0],[0,-1]]);
one:=Matrix([[1,0],[0,1]]);
zero:=Matrix([[0,0],[0,0]]);

make_spin_matrix:=proc(blockmat)
 local UL,UR,LL,LR:
 UL:=blockmat[1][1]:
 UR:=blockmat[1][2]:
 LL:=blockmat[2][1]:
 LR:=blockmat[2][2]:
 return Matrix([[UL[1,1],UL[1,2],UR[1,1],UR[1,2]],
                [UL[2,1],UL[2,2],UR[2,1],UR[2,2]],
                [LL[1,1],LL[1,2],LR[1,1],LR[1,2]],
                [LL[2,1],LL[2,2],LR[2,1],LR[2,2]]]):
end:

is_zero:=proc(mat)
 local flag,i,j:
 flag:=true:
 for i from 1 to 4 do
 for j from 1 to 4 do
    flag:=flag and evalb((mat[i,j]=0)):
    end do:
    end do:
 return flag:
end:


check_gamma_matrices:=proc(g1,g2,g3,g4,g5)
 local flag:
 flag:=is_zero(g1.g2+g2.g1):
 flag:=flag and is_zero(g1.g3+g3.g1):
 flag:=flag and is_zero(g1.g4+g4.g1):
 flag:=flag and is_zero(g2.g3+g3.g2):
 flag:=flag and is_zero(g2.g4+g4.g2):
 flag:=flag and is_zero(g3.g4+g4.g3):
 flag:=flag and is_zero(g1.g1-IdentityMatrix(4,4)):
 flag:=flag and is_zero(g2.g2-IdentityMatrix(4,4)):
 flag:=flag and is_zero(g3.g3-IdentityMatrix(4,4)):
 flag:=flag and is_zero(g4.g4-IdentityMatrix(4,4)):
 flag:=flag and is_zero(g5.g5-IdentityMatrix(4,4)):
 flag:=flag and is_zero(g5-g4.g1.g2.g3):
 return flag:
end:


gamGR1:=make_spin_matrix([[zero,I*sigma1],[-I*sigma1,zero]]);
gamGR2:=make_spin_matrix([[zero,-I*sigma2],[I*sigma2,zero]]);
gamGR3:=make_spin_matrix([[zero,I*sigma3],[-I*sigma3,zero]]);
gamGR4:=make_spin_matrix([[zero,one],[one,zero]]);
gamGR5:=make_spin_matrix([[-one,zero],[zero,one]]);
print("Checking gamGR matrices");
check_gamma_matrices(gamGR1,gamGR2,gamGR3,gamGR4,gamGR5);

#  conversion matrices to go to Degrand-Rossi

C_DP_to_GR:=I/sqrt(2)*make_spin_matrix([[sigma2,-sigma2],[sigma2,sigma2]]);
C_UK_to_GR:=I/sqrt(2)*make_spin_matrix([[-sigma2,-sigma2],[-sigma2,sigma2]]);
C_CH_to_GR:=I*make_spin_matrix([[0,-sigma2],[sigma2,0]]);

#  check Dirac-Pauli gamma matrices
# gam1=((0,-i*s1)(i*s1,0)) gam2=((0,-i*s2)(i*s2,0)) gam3=((0,-i*s3)(i*s3,0)) gam4=((1,0)(0,-1)) gam5=((0,1)(1,0))

gamDP1:=make_spin_matrix([[zero,-I*sigma1],[I*sigma1,zero]]); 
gamDP2:=make_spin_matrix([[zero,-I*sigma2],[I*sigma2,zero]]); 
gamDP3:=make_spin_matrix([[zero,-I*sigma3],[I*sigma3,zero]]); 
gamDP4:=make_spin_matrix([[one,zero],[zero,-one]]); 
gamDP5:=make_spin_matrix([[zero,one],[one,zero]]);
print("Checking gamDP matrices");
check_gamma_matrices(gamDP1,gamDP2,gamDP3,gamDP4,gamDP5);

is_zero(gamDP1-MatrixInverse(C_DP_to_GR).gamGR1.C_DP_to_GR);
is_zero(gamDP2-MatrixInverse(C_DP_to_GR).gamGR2.C_DP_to_GR);
is_zero(gamDP3-MatrixInverse(C_DP_to_GR).gamGR3.C_DP_to_GR);
is_zero(gamDP4-MatrixInverse(C_DP_to_GR).gamGR4.C_DP_to_GR);
is_zero(gamDP5-MatrixInverse(C_DP_to_GR).gamGR5.C_DP_to_GR);

#  check UKQCD gamma matrices
# gam1=((0,i*s1)(-i*s1,0)) gam2=((0,i*s2)(-i*s2,0)) gam3=((0,i*s3)(-i*s3,0)) gam4=((1,0)(0,-1)) gam5=((0,-1)(-1,0))

gamUK1:=make_spin_matrix([[zero,I*sigma1],[-I*sigma1,zero]]); 
gamUK2:=make_spin_matrix([[zero,I*sigma2],[-I*sigma2,zero]]); 
gamUK3:=make_spin_matrix([[zero,I*sigma3],[-I*sigma3,zero]]); 
gamUK4:=make_spin_matrix([[one,zero],[zero,-one]]); 
gamUK5:=make_spin_matrix([[zero,-one],[-one,zero]]);
print("Checking gamUK matrices");
check_gamma_matrices(gamUK1,gamUK2,gamUK3,gamUK4,gamUK5);

is_zero(gamUK1-MatrixInverse(C_UK_to_GR).gamGR1.C_UK_to_GR);
is_zero(gamUK2-MatrixInverse(C_UK_to_GR).gamGR2.C_UK_to_GR);
is_zero(gamUK3-MatrixInverse(C_UK_to_GR).gamGR3.C_UK_to_GR);
is_zero(gamUK4-MatrixInverse(C_UK_to_GR).gamGR4.C_UK_to_GR);
is_zero(gamUK5-MatrixInverse(C_UK_to_GR).gamGR5.C_UK_to_GR);

#  check chiral gamma matrices
# gam1=((0,-i*s1)(i*s1,0)) gam2=((0,-i*s2)(i*s2,0)) gam3=((0,-i*s3)(i*s3,0)) gam4=((0,-1)(-1,0))  gam5=((1,0)(0,-1))

gamCH1:=make_spin_matrix([[zero,-I*sigma1],[I*sigma1,zero]]); 
gamCH2:=make_spin_matrix([[zero,-I*sigma2],[I*sigma2,zero]]); 
gamCH3:=make_spin_matrix([[zero,-I*sigma3],[I*sigma3,zero]]); 
gamCH4:=make_spin_matrix([[zero,-one],[-one,zero]]);  
gamCH5:=make_spin_matrix([[one,zero],[zero,-one]]);
print("Checking gamCH matrices");
check_gamma_matrices(gamCH1,gamCH2,gamCH3,gamCH4,gamCH5);
 
is_zero(gamCH1-MatrixInverse(C_CH_to_GR).gamGR1.C_CH_to_GR);
is_zero(gamCH2-MatrixInverse(C_CH_to_GR).gamGR2.C_CH_to_GR);
is_zero(gamCH3-MatrixInverse(C_CH_to_GR).gamGR3.C_CH_to_GR);
is_zero(gamCH4-MatrixInverse(C_CH_to_GR).gamGR4.C_CH_to_GR);
is_zero(gamCH5-MatrixInverse(C_CH_to_GR).gamGR5.C_CH_to_GR);


C_CH_to_UK:=(1/sqrt(2))*make_spin_matrix([[-one,one],[one,one]]);
is_zero(gamCH1-MatrixInverse(C_CH_to_UK).gamUK1.C_CH_to_UK);
is_zero(gamCH2-MatrixInverse(C_CH_to_UK).gamUK2.C_CH_to_UK);
is_zero(gamCH3-MatrixInverse(C_CH_to_UK).gamUK3.C_CH_to_UK);
is_zero(gamCH4-MatrixInverse(C_CH_to_UK).gamUK4.C_CH_to_UK);
is_zero(gamCH5-MatrixInverse(C_CH_to_UK).gamUK5.C_CH_to_UK);

MatrixInverse(C_UK_to_GR).C_DP_to_GR;


Ferm:=Vector([a[0],a[1],a[2],a[3]]);
