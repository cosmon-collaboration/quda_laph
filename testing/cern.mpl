interface(quiet=true):

Nx:=12:  Ny:=12:  Nz:=12:  Nt:=16:
Px:=3:   Py:=4:   Pz:=2:   Pt:=4:

Rx:=Nx/Px: Ry:=Ny/Py: Rz:=Nz/Pz: Rt:=Nt/Pt:

get_comm_coord:=proc(x,y,z,t)
 global Nx,Ny,Nz,Nt,Rx,Ry,Rz,Rt:
 local Cx,Cy,Cz,Ct:
 return [iquo(x,Rx),iquo(y,Ry),iquo(z,Ry),iquo(t,Rt)]:
end: 

for ccx from 0 to Px-1 do
for ccy from 0 to Py-1 do
for ccz from 0 to Pz-1 do
for cct from 0 to Pt-1 do
   count:=0:
   lexico:=0:
   printf("\n\nComm coordinate [%2d,%2d,%2d,%2d]\n\n",ccx,ccy,ccz,cct):
   for t from 0 to Nt-1 do
   for x from 0 to Nx-1 do
   for y from 0 to Ny-1 do
   for z from 0 to Nz-1 do
      if type((x+y+z+t),odd) then
         cc:=get_comm_coord(x,y,z,t):
         if (cc[1]=ccx)and(cc[2]=ccy)and(cc[3]=ccz)and(cc[4]=cct) then
            printf("Count: %4d Lexico: %4d   Coordinate: [%2d,%2d,%2d,%2d]\n",count,lexico,x,y,z,t):
            count:=count+1:
         end if:
         lexico:=lexico+1:
      end if:
      end do:
      end do:
      end do:
      end do:
   end do:
   end do:
   end do:
   end do:
