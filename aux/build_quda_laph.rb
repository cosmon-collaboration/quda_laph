#!/usr/bin/ruby

require "pathname.rb"

$: << File.expand_path(File.dirname(Pathname.new(__FILE__).realpath))

require "fileutils.rb"
require "date.rb"


class Launcher

 #############################################
 #
 #       initialize Launcher
 #
 #############################################

 def initialize

   # determine the machine type

  if File.exists?("machine_name")
     @machine=File.read("machine_name")
  else
     @machine=`hostname -f`
  end

  if (@machine.index("qcdcomp-2-1")!=nil)
     @machine="cmucluster"
  elsif (@machine.index("bridges")!=nil)
     @machine="bridges"
  elsif (@machine.index("perlmutter")!=nil)
     @machine="perlmutter"
  elsif (@machine.index("summit")!=nil)
     @machine="summit"
  elsif (@machine.index("frontier")!=nil)
     @machine="frontier"
  else
     abort(" Unsupported machine "+`hostname`)
  end

  if (@machine=="bridges")

     puts
     puts "  intelmpi module must be used"
     puts "    must load modules  phdf5/hdf5"
     puts "  Are these module loaded? (y/n)"
     reply=gets.chomp.strip
     abort("aborting...") if (reply!="y")
     puts

  end

  if (@machine=="perlmutter")

     puts " Will use cray compiler; are modules below set? "
     puts "   module load cpu"
     puts "   module load PrgEnv-cray"
     puts "   module unload darshan"
     puts "   module load cray-hdf5-parallel"
     puts "  Are these module loaded? (y/n)"

  end

  if (@machine=="summit")

     puts
     puts "  IBM compiler fails/takes very long to compile chroma"
     puts "  gcc/12.1.0, hdf5, spectrum-mpi, netlib-lapack, libxml2, binutils, gmp modules must be loaded"
     puts "  Use an interactive bsub to compile on a node: then can use -j 32"
     puts "  Are these modules loaded? (y/n)"
     reply=gets.chomp.strip
     abort("aborting...") if (reply!="y")
     puts

  end

  puts
  puts " QUDA_LAPH software installation on machine '"+@machine+"'"
  puts
  
      # default directories

  now=DateTime.now
  nowstr="DATE-"+now.year.to_s+"-"+now.month.to_s+"-"+now.day.to_s

  if (@machine=="cmucluster")
     @qudalaphdir="/home/colin/quda_laph"
     @srcdir=@qudalaphdir+"/source"
     @builddir=@qudalaphdir+"/build/"+nowstr
     @installdir=@qudalaphdir+"/install/"+nowstr
     @qudainstalldir="/home/colin/quda/install/DATE-2024-10-30"
  elsif (@machine=="perlmutter")
     @qudalaphdir="/global/cfs/cdirs/m2986/stoch_nn/colin/quda_laph"
     @srcdir=@qudalaphdir+"/source"
     @builddir=@qudalaphdir+"/build/"+nowstr
     @installdir=@qudalaphdir+"/install/"+nowstr
     @qudainstalldir="/global/cfs/cdirs/m2986/stoch_nn/colin/quda/install/DATE-2024-10-18"
  elsif (@machine=="summit")
     @qudalaphdir="/autofs/nccs-svm1_proj/nph162/colin/quda_laph"
     @srcdir=@qudalaphdir+"/source"
     @builddir=@qudalaphdir+"/build/"+nowstr
     @installdir=@qudalaphdir+"/install/"+nowstr
     @qudainstalldir=""
  elsif (@machine=="frontier")
     @qudalaphdir="/autofs/nccs-svm1_proj/nph165/colin/quda_laph"
     @srcdir=@qudalaphdir+"/source"
     @builddir=@qudalaphdir+"/build/"+nowstr
     @installdir=@qudalaphdir+"/install/"+nowstr
     @qudainstalldir="/autofs/nccs-svm1_proj/nph162/colin/quda/install/DATE-2024-10-18"
  elsif (@machine=="bridges")
     @qudalaphdir="/jet/home/mornings/quda_laph"
     @srcdir=@qudalaphdir+"/source"
     @builddir=@qudalaphdir+"/build/"+nowstr
     @installdir=@qudalaphdir+"/install/"+nowstr
     @qudainstalldir=""
  else
     abort("Unsupported machine ")
  end

      # set the QUDA_LAPH source directory

  puts " The default QUDA_LAPH <srcdir> directory for this machine is "
  puts
  puts "        "+@srcdir
  puts
  puts " Hit <enter> to accept the default, or enter a different path: "
  print " QUDA_LAPH source directory? "
  reply=gets.chomp.strip
  if (reply!="")
     @srcdir=reply
     puts " The current QUDA_LAPH source directory is "
     puts "        "+@srcdir
  end
  
      # set architecture mpi or serial

  puts
  puts
  commchoices=Array.new
  commchoices.push("mpi-parallel")
  commchoices.push("none-serial")
  puts "Which type of communication to use?"
  nchoices=commchoices.size
  for k in 1..nchoices
     puts " ["+k.to_s+"]  "+commchoices[k-1]
     end
  print "Choice? "
  choice = gets.to_i
  abort("Invalid entry") if (choice<1)||(choice>nchoices)
  @comm=commchoices[choice-1]
  commstr = { "mpi-parallel" => "mpi", "none-serial" => "ser" }
  suffix = "_"+commstr[@comm]

  @builddir=@builddir+suffix
  @installdir=@installdir+suffix
  @qudainstalldir=@qudainstalldir+suffix
  if (@comm=="mpi-parallel")
     @arch=" -DARCH=PARALLEL"
  else
     @arch=" -DARCH=SERIAL"
  end

      # set the build directory

  puts
  puts " The default build directory for this machine is "
  puts
  puts "        "+@builddir
  puts
  puts " Hit <enter> to accept the default, or enter a different path: "
  print " Build directory? "
  reply=gets.chomp.strip
  if (reply!="")
     @builddir=reply
     puts " The current build directory is "
     puts "        "+@builddir
  end

      # set the installation directory

  puts
  puts " The default installation directory for this machine is "
  puts
  puts "        "+@installdir
  puts
  puts " Hit <enter> to accept the default, or enter a different path: "
  print " Installation directory? "
  reply=gets.chomp.strip
  if (reply!="")
     @installdir=reply
     puts " The current installation directory is "
     puts "        "+@installdir
  end

      # set the QUDA installation directory

  puts
  puts " The default QUDA installation directory for this machine is "
  puts
  puts "        "+@qudainstalldir
  puts
  puts " Hit <enter> to accept the default, or enter a different path: "
  print " QUDA installation directory? "
  reply=gets.chomp.strip
  if (reply!="")
     @qudainstalldir=reply
     puts " The current QUDA installation directory is "
     puts "        "+@qudainstalldir
  end

 end


 def checksource

     # check for existence of source code directories

  abort("QUDA_LAPH source directory does not exist") if !((File.exist?(@srcdir))&&(File.directory?(@srcdir)))    
  
 end


 def setup

  @makeopt=""
 
  if (@machine=="cmucluster")

     if (@comm=="mpi-parallel")
        @cc="/usr/lib64/openmpi/bin/mpicc"
        @cxx="/usr/lib64/openmpi/bin/mpicxx"
     elsif (@comm=="none-serial")
        @cc="/usr/bin/gcc"
        @cxx="/usr/bin/g++"
     end
     @gpu="cuda"
     @gpudir="/usr/local/cuda-11.4"
     #@gpubindir=@gpudir+"/bin"
     #@gpuincdir=@gpudir+"/include"
     #@gpuarch="sm_86"
     @makeopt="-j16"
     @cflags="-O3 -fopenmp"
     @cxxflags="-O3 -fopenmp"
     @linkerflags=""
     @hostcblas="gsl"
     @cblaslib="/usr/lib64/libgslcblas.so"
     @cblasinc="/usr/include/gsl"
#     @hostcblas="openblas"
#     @cblaslib="/usr/lib64/libopenblas64.so"
#     @cblasinc="/usr/include/openblas"

  elsif (@machine=="perlmutter")
  
     if (@comm=="mpi-parallel")
        @cc="cc"
        @cxx="CC"
     elsif (@comm=="none-serial")
        abort("Not currently supported")
     end
     @gpu="cuda"
     @gpudir="/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/cuda/12.2"
     #@gpubindir=@gpudir+"/bin"
     #@gpuincdir=@gpudir+"/include"
     #@gpuarch="sm_80"
     @makeopt="-j4"
     @cflags="-O3 -fopenmp"
     @cxxflags="-O3 -fopenmp"
     @linkerflags=""
     @hostcblas="gsl"
     @cblaslib=""
     @cblasinc=""
          
  elsif (@machine=="bridges")
  
     if (@comm=="mpi-parallel")
        @cc="mpiicc"
        @cxx="mpiicpc"
     elsif (@comm=="none-serial")
        @cc="icc"
        @cxx="icpc"
     end
     @gpu="cuda"
     @gpudir="/opt/packages/cuda/v12.4.0"
     #@gpubindir=@gpudir+"/bin"
     #@gpuincdir=@gpudir+"/include"
     #@gpuarch="sm_70"
     @makeopt="-j16"
     @cflags="-O3 -qopenmp"
     @cxxflags="-O3 -qopenmp"
     @linkerflags=""
     @hostcblas="gsl"
     @cblaslib=""
     @cblasinc=""

  elsif (@machine=="summit")
  
     if (@comm=="mpi-parallel")||(@comm=="qmp-parallel")
        @cc="mpicc"
        @cxx="mpiCC"
     elsif (@comm=="none-serial")
        abort("Not currently supported")
     end
     @gpu="cuda"
     @gpudir="/sw/summit/cuda/12.2.0"
     #@gpubindir=@gpudir+"/bin"
     #@gpuincdir=@gpudir+"/include"
     #@gpuarch="sm_70"
     @makeopt="-j32"
     @cflags="-O3 -fopenmp"
     @cxxflags="-O3 -fopenmp"
     @linkerflags=""
     @hostcblas="gsl"
     @cblaslib=""
     @cblasinc=""

  elsif (@machine=="frontier")
  
     if (@comm=="mpi-parallel")||(@comm=="qmp-parallel")
        @cc="cc"
        @cxx="CC"
     elsif (@comm=="none-serial")
        abort("Not currently supported")
     end
     @gpu="hip"
     @gpudir="/opt/rocm-5.7.1"
     #@gpubindir=@gpudir+"/bin"
     #@gpuincdir=@gpudir+"/include"
     #@gpuarch="gfx90a"
     @makeopt="-j32"
     @cflags=" -D__HIP_PLATFORM_AMD__ -O2  -fopenmp"
     @cxxflags=" -D__HIP_PLATFORM_AMD__ -O2 -fopenmp"
     @linkerflags=""
     @hostcblas="openblas"
     @cblaslib="/sw/frontier/spack-envs/core-24.07/opt/gcc-7.5.0/openblas-0.3.26-3m42udfwiyrshafp5qag4e2fixa4pnbq/lib/libopenblas.so"
     @cblasinc="/sw/frontier/spack-envs/core-24.07/opt/gcc-7.5.0/openblas-0.3.26-3m42udfwiyrshafp5qag4e2fixa4pnbq/include"

  end
 end


 def buildquda

   puts
   puts "Build QUDA_LAPH? (y/n) "
   reply=gets.chomp.strip
   return if (reply!="y")
   
   FileUtils.rm_rf(@installdir)
   FileUtils.rm_rf(@builddir)
   FileUtils.mkdir_p @builddir

   if (@hostcblas=="gsl")
      hostcblas="GSL"
   elsif (@hostcblas=="openblas")
      hostcblas="OPENBLAS"
   else
      hostcblas="NONE"
   end

#  
#cmd="cmake --fresh  -S .. -B . -DARCH=PARALLEL"+
#   " -DCMAKE_CXX_COMPILER=\"CC\""+
#   " -DCMAKE_C_COMPILER=\"cc\""+
#   " -DGPUTOOLKIT_DIR=/opt/rocm-5.3.0"+
#   " -DQUDA_DIR=/autofs/nccs-svm1_proj/nph162/colin/quda/install/DATE-2024-10-18_mpi"+
#   " -DBUILD_TESTS=ON"+
#   " -DHOSTCBLAS=OPENBLAS"+
#   " -DCBLAS_INC=/sw/frontier/spack-envs/core-24.07/opt/gcc-7.5.0/openblas-0.3.26-3m42udfwiyrshafp5qag4e2fixa4pnbq/include"+
#   " -DCBLAS_LIB=/sw/frontier/spack-envs/core-24.07/opt/gcc-7.5.0/openblas-0.3.26-3m42udfwiyrshafp5qag4e2fixa4pnbq/lib/libopenblas.so"+
#   " -DCMAKE_CXX_STANDARD=17"+
#   " -DCMAKE_CXX_FLAGS=\" -D__HIP_PLATFORM_AMD__ -O2 -Wall -fopenmp\""+
#   " -DCMAKE_C_FLAGS=\" -D__HIP_PLATFORM_AMD__ -O2 -Wall -fopenmp\""



  cmd="cmake --fresh -S "+@srcdir+" -B "+@builddir+" "+@arch+
      " -DCMAKE_CXX_COMPILER="+@cxx+
      " -DCMAKE_C_COMPILER="+@cc+
      " -DCMAKE_EXE_LINKER_FLAGS=\""+@linkerflags+"\""+
      " -DCMAKE_INSTALL_PREFIX="+@installdir+
      " -DCMAKE_CXX_STANDARD=17"+
      " -DCMAKE_CXX_EXTENSIONS=OFF"+
      " -DCMAKE_CXX_FLAGS=\""+@cxxflags+"\""+
      " -DCMAKE_C_FLAGS=\""+@cflags+"\""+
      " -DQUDA_DIR=\""+@qudainstalldir+"\""+
      " -DGPUTOOLKIT_DIR=\""+@gpudir+"\""+
      " -DHOSTCBLAS=\""+hostcblas+"\""+
      " -DCBLAS_INC=\""+@cblasinc+"\""+
      " -DCBLAS_LIB=\""+@cblaslib+"\""



#      gpuoptions+
#      " -DQUDA_DIRAC_CLOVER=ON"+
#      " -DQUDA_DIRAC_DOMAIN_WALL=OFF"+
#      " -DQUDA_MDW_FUSED_LS_LIST=\"8,12,16,20\""+
#      " -DQUDA_DIRAC_STAGGERED=OFF"+
#      " -DQUDA_DIRAC_TWISTED_MASS=OFF"+
#      " -DQUDA_DIRAC_TWISTED_CLOVER=OFF"+
#      " -DQUDA_DIRAC_WILSON=ON"+
#      " -DQUDA_FORCE_GAUGE=OFF"+
#      " -DQUDA_FORCE_HISQ=OFF"+
#      " -DQUDA_GAUGE_ALG=OFF"+
#      " -DQUDA_GAUGE_TOOLS=OFF"+
#      " -DQUDA_CLOVER_DYNAMIC=ON"+
#      " -DQUDA_CLOVER_RECONSTRUCT=ON"+
#      " -DQUDA_QDPJIT=OFF"+
#      " -DQUDA_INTERFACE_QDPJIT=OFF"+
#      " -DQUDA_INTERFACE_MILC=OFF"+
#      " -DQUDA_INTERFACE_CPS=OFF"+
#      " -DQUDA_INTERFACE_QDP=ON"+
#      " -DQUDA_INTERFACE_TIFR=OFF"+
#      " -DQUDA_QMP="+qmp+
#      " -DQUDA_QIO=OFF"+
#      " -DQUDA_OPENMP=ON"+
#      " -DQUDA_MULTIGRID=ON"+
#      " -DQUDA_LAPLACE=ON"+
#      " -DQUDA_MPI="+mpi+
#      " -DQUDA_NVSHMEM=OFF"+
#      " -DQUDA_MAX_MULTI_BLAS_N=9"+
#      " -DQUDA_DOWNLOAD_EIGEN=ON"+
#      " -DQUDA_EIGEN_VERSION=3.4.0"+
#      " -DQUDA_DOWNLOAD_USQCD=OFF"+
#      " -DCMAKE_BUILD_TYPE=""RELEASE"""+
#      " -DBUILD_SHARED_LIBS=ON"+
#      " -DQUDA_BUILD_SHAREDLIB=ON"+
#      " -DQUDA_BUILD_ALL_TESTS=ON"+
#      " -DQUDA_CTEST_DISABLE_BENCHMARKS=OFF"+
#      " -DQUDA_BACKWARDS=ON"

   puts cmd

     # do cmake
   system(cmd)
   success=($?.to_i==0)
   abort("QUDA_LAPH cmake failed; aborting") if (!success)
     # do make
   system("cmake --build "+@builddir+" "+@makeopt+" 2>&1 | tee "+@builddir+"/build.log")
   success=($?.to_i==0)
   abort("QUDA_LAPH build failed; aborting") if (!success)

 end
 

 def launch
   checksource
   setup
   buildquda
 end
    
end


#####################################################################
#
#                           MAIN DRIVER
#
#####################################################################

main = Launcher.new

main.launch

