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
  puts " QUDA software installation on machine '"+@machine+"'"
  puts
  
      # default directories

  now=DateTime.now
  nowstr="DATE-"+now.year.to_s+"-"+now.month.to_s+"-"+now.day.to_s

  if (@machine=="cmucluster")
     @qudadir="/home/colin/quda"
     @srcdir=@qudadir+"/source"
     @builddir=@qudadir+"/build/"+nowstr
     @installdir=@qudadir+"/install/"+nowstr
  elsif (@machine=="perlmutter")
     @qudadir="/global/cfs/cdirs/m2986/stoch_nn/colin/quda"
     @srcdir="/global/cfs/cdirs/m2986/stoch_nn/colin/quda/source"
     @builddir=@qudadir+"/build/"+nowstr
     @installdir="/global/cfs/cdirs/m2986/stoch_nn/colin/quda/install/perlmutter/"+nowstr
  elsif (@machine=="summit")
     @qudadir="/autofs/nccs-svm1_proj/nph162/colin/quda"
     @srcdir="/autofs/nccs-svm1_proj/nph162/colin/quda/source"
     @builddir=@qudadir+"/build/"+nowstr
     @installdir="/autofs/nccs-svm1_proj/nph162/colin/quda/install/"+nowstr
  elsif (@machine=="frontier")
     @qudadir="/autofs/nccs-svm1_proj/nph165/colin/quda"
     @srcdir="/autofs/nccs-svm1_proj/nph165/colin/quda/source"
     @builddir=@qudadir+"/build/"+nowstr
     @installdir="/autofs/nccs-svm1_proj/nph165/colin/quda/install/"+nowstr
  elsif (@machine=="bridges")
     @qudadir="/jet/home/mornings/quda"
     @srcdir=@qudadir+"/source"
     @builddir=@qudadir+"/build/"+nowstr
     @installdir=@qudadir+"/install/"+nowstr
  else
     abort("Unsupported machine ")
  end

      # set the QUDA source directory

  puts " The directory structure for the QUDA source code should be:"
  puts " The default QUDA <srcdir> directory for this machine is "
  puts
  puts "        "+@srcdir
  puts
  puts " Hit <enter> to accept the default, or enter a different path: "
  print " QUDA source directory? "
  reply=gets.chomp.strip
  if (reply!="")
     @srcdir=reply
     puts " The current QUDA source directory is "
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

 end

 
 def checksource

     # check for existence of source code directories

  abort("QUDA source directory does not exist") if !((File.exist?(@srcdir))&&(File.directory?(@srcdir)))    
  
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
     @gpubindir=@gpudir+"/bin"
     @gpuincdir=@gpudir+"/include"
     @gpuarch="sm_86"
     @makeopt="-j8"
     @cflags="-O3 -fopenmp"
     @cxxflags="-O3 -fopenmp"
     @linkerflags=""

  elsif (@machine=="perlmutter")
  
     if (@comm=="mpi-parallel")
        @cc="cc"
        @cxx="CC"
     elsif (@comm=="none-serial")
        abort("Not currently supported")
     end
     @gpu="cuda"
     @gpudir="/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/cuda/12.2"
     @gpubindir=@gpudir+"/bin"
     @gpuincdir=@gpudir+"/include"
     @gpuarch="sm_80"
     @makeopt="-j4"
     @cflags="-O3 -fopenmp"
     @cxxflags="-O3 -fopenmp"
     @linkerflags=""
          
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
     @gpubindir=@gpudir+"/bin"
     @gpuincdir=@gpudir+"/include"
     @gpuarch="sm_70"
     @makeopt="-j8"
     @cflags="-O3 -qopenmp"
     @cxxflags="-O3 -qopenmp"
     @linkerflags=""

  elsif (@machine=="summit")
  
     if (@comm=="mpi-parallel")
        @cc="mpicc"
        @cxx="mpiCC"
     elsif (@comm=="none-serial")
        abort("Not currently supported")
     end
     @gpu="cuda"
     @gpudir="/sw/summit/cuda/12.2.0"
     @gpubindir=@gpudir+"/bin"
     @gpuincdir=@gpudir+"/include"
     @gpuarch="sm_70"
     @makeopt="-j8"
     @cflags="-O3 -fopenmp"
     @cxxflags="-O3 -fopenmp"
     @linkerflags=" -lsframe"

  elsif (@machine=="frontier")
  
     if (@comm=="mpi-parallel")
        @cc="hipcc"
        @cxx="hipcc"
     elsif (@comm=="none-serial")
        abort("Not currently supported")
     end
     @gpu="hip"
     @gpudir="/opt/rocm-6.0.0"
     @gpubindir=@gpudir+"/bin"
     @gpuincdir=@gpudir+"/include"
     @gpuarch="gfx90a"
     @makeopt="-j8"
     @cflags="-O3 -fopenmp  -D__HIP_PLATFORM_AMD__ -Wall --offload-arch=gfx90a"
     @cxxflags="-O3 -fopenmp  -D__HIP_PLATFORM_AMD__ -Wall --offload-arch=gfx90a"
     @linkerflags="--offload-arch=gfx90a"

  end
 end


 def buildquda

   puts
   puts "Build QUDA? (y/n) "
   reply=gets.chomp.strip
   return if (reply!="y")
   
   FileUtils.rm_rf(@installdir)
   FileUtils.rm_rf(@builddir)
   FileUtils.mkdir_p @builddir

   if (@comm=="mpi-parallel")
      mpi="ON"
   else
      mpi="OFF"
   end

   if (@gpu=="cuda")
      gpuoptions=" -DQUDA_TARGET_TYPE=\"CUDA\""+
                 " -DQUDA_GPU_ARCH="+@gpuarch+
                 " -DCUDAToolkit_BIN_DIR="+@gpubindir+
                 " -DCUDAToolkit_CUPTI_INCLUDE_DIR="+@gpuincdir
   elsif (@gpu=="hip")
      gpuoptions=" -DQUDA_TARGET_TYPE=\"HIP\""+
                 " -DQUDA_GPU_ARCH="+@gpuarch+
                 " -DROCM_PATH=\""+@gpudir+"\""
   elsif (@gpu=="sysl")
      abort("Sysl not yet supported")
   end
   
   cmd="cmake -S "+@srcdir+" -B "+@builddir+
       " -DCMAKE_CXX_COMPILER="+@cxx+
       " -DCMAKE_C_COMPILER="+@cc+
       " -DCMAKE_INSTALL_PREFIX="+@installdir+
       " -DCMAKE_CXX_STANDARD=17"+
       " -DCMAKE_CXX_EXTENSIONS=OFF"+
       " -DCMAKE_CXX_FLAGS=\""+@cxxflags+"\""+
       " -DCMAKE_C_FLAGS=\""+@cflags+"\""+
       gpuoptions+
       " -DQUDA_DIRAC_CLOVER=ON"+
       " -DQUDA_DIRAC_DOMAIN_WALL=OFF"+
       " -DQUDA_DIRAC_CLOVER_HASENBUSCH=OFF"+
       " -DQUDA_MDW_FUSED_LS_LIST=\"8,12,16,20\""+
       " -DQUDA_DIRAC_STAGGERED=OFF"+
       " -DQUDA_DIRAC_TWISTED_MASS=OFF"+
       " -DQUDA_DIRAC_TWISTED_CLOVER=OFF"+
       " -DQUDA_DIRAC_WILSON=ON"+
       " -DQUDA_FORCE_GAUGE=OFF"+
       " -DQUDA_FORCE_HISQ=OFF"+
       " -DQUDA_GAUGE_ALG=OFF"+
       " -DQUDA_GAUGE_TOOLS=OFF"+
       " -DQUDA_CLOVER_DYNAMIC=ON"+
       " -DQUDA_CLOVER_RECONSTRUCT=ON"+
       " -DQUDA_QDPJIT=OFF"+
       " -DQUDA_INTERFACE_QDPJIT=OFF"+
       " -DQUDA_INTERFACE_MILC=OFF"+
       " -DQUDA_INTERFACE_CPS=OFF"+
       " -DQUDA_INTERFACE_QDP=ON"+
       " -DQUDA_INTERFACE_TIFR=OFF"+
       " -DQUDA_QMP=OFF"+
       " -DQUDA_QIO=OFF"+
       " -DQUDA_OPENMP=ON"+
       " -DQUDA_MULTIGRID=ON"+
       " -DQUDA_LAPLACE=ON"+
       " -DQUDA_MPI="+mpi+
       " -DQUDA_NVSHMEM=OFF"+
       " -DQUDA_MAX_MULTI_BLAS_N=9"+
       " -DQUDA_DOWNLOAD_EIGEN=ON"+
       " -DQUDA_EIGEN_VERSION=3.4.0"+
       " -DQUDA_DOWNLOAD_USQCD=OFF"+
       " -DCMAKE_BUILD_TYPE=""RELEASE"""+
       " -DBUILD_SHARED_LIBS=ON"+
       " -DQUDA_BUILD_SHAREDLIB=ON"+
       " -DQUDA_BUILD_ALL_TESTS=ON"+
       " -DQUDA_CTEST_DISABLE_BENCHMARKS=OFF"+
       " -DQUDA_BACKWARDS=ON"+
       " -DCMAKE_EXE_LINKER_FLAGS=\""+@linkerflags+"\""

   puts cmd

     # do cmake
   system(cmd)
   success=($?.to_i==0)
   abort("QUDA cmake failed; aborting") if (!success)
     # do make
   system("cmake --build "+@builddir+" "+@makeopt+" 2>&1 | tee "+@builddir+"/build.log")
   success=($?.to_i==0)
   abort("QUDA build failed; aborting") if (!success)
     # do make install
   system("cmake --install "+@builddir)
   success=($?.to_i==0)
   abort("QUDA build failed; aborting") if (!success)

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

