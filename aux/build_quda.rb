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
  elsif (@machine.index("qcdcomp-2-2")!=nil)
     @machine="cmucluster"
  elsif (@machine.index("bridges")!=nil)
     @machine="bridges"
  elsif (@machine.index("perlmutter")!=nil)
     @machine="perlmutter"
  elsif (@machine.index("summit")!=nil)
     @machine="summit"
  elsif (@machine.index("frontier")!=nil)
     @machine="frontier"
  elsif (@machine.index("tioga")!=nil)
     @machine="tioga"
  else
     abort(" Unsupported machine "+`hostname`)
  end

  puts
  puts " QUDA software installation on machine '"+@machine+"'"
  puts
  cmds=[]
  vsuffix = ""

  if (@machine=="bridges")

     puts " Do the modules need to be set? (y/n) "
     reply=gets.chomp.strip
     if (reply=="y")
        cmds=["module purge",
              "module load allocations",
              "module load gcc/10.2.0",
              "module load openmpi/4.0.5-gcc10.2.0",
              "module load phdf5/1.10.7-openmpi4.0.5-gcc10.2.0",
              "module load cuda/12.6.1",
              "module load openblas/0.3.12-gcc10.2.0"]
     end
  end

  if (@machine=="perlmutter")

     puts " Do the modules need to be set? (y/n) "
     reply=gets.chomp.strip
     if (reply=="y")
        cmds=["module purge",
              "module load PrgEnv-gnu/8.5.0",
              "module load gpu",
              "module load cmake/3.30.2",
              "module load cpe-cuda/23.12",
              "module load cudatoolkit/12.2",
              "module load craype-accel-nvidia80",
              "module load craype-x86-milan",
              "module load eigen/3.4.0"]
     end
  end

  if (@machine=="summit")

     puts " Do the modules need to be set? (y/n) "
     reply=gets.chomp.strip
     if (reply=="y")
        cmds=["module purge",
              "module load gcc/12.1.0",
              "module load cuda/12.2.0",
              "module load spectrum-mpi/10.4.0.6-20230210",
              "module load binutils/2.40",
              "module load eigen/3.4.0"]
     end
  end

  if (@machine=="frontier")

     envchoices=Array.new
     envchoices.push("cray")
     envchoices.push("amd")
     puts "Which environment to use?"
     nchoices=envchoices.size
     for k in 1..nchoices
        puts " ["+k.to_s+"]  "+envchoices[k-1]
        end
     print "Choice? "
     choice = gets.to_i
     abort("Invalid entry") if (choice<1)||(choice>nchoices)
     @env=envchoices[choice-1]
     envstr = { "cray" => "cray", "amd" => "amd" }
     vsuffix += "_"+envstr[@env]

     rocmchoices=Array.new
     rocmchoices.push("5.3.0")
     rocmchoices.push("5.4.0")
     rocmchoices.push("5.4.3") 
     rocmchoices.push("5.5.1")
     rocmchoices.push("5.6.0")
     rocmchoices.push("5.7.0")
     rocmchoices.push("5.7.1")
     rocmchoices.push("6.2.4")
     rocmchoices.push("6.3.1")
     puts "Which rocm version?"
     nchoices=rocmchoices.size
     for k in 1..nchoices
        puts " ["+k.to_s+"]  "+rocmchoices[k-1]
        end
     print "Choice? "
     choice = gets.to_i
     abort("Invalid entry") if (choice<1)||(choice>nchoices)
     @rocm=rocmchoices[choice-1]
     rocmstr = { "5.3.0" => "530", "5.4.0" => "540", "5.4.3" => "543", 
                 "5.5.1" => "551", "5.6.0" => "560", "5.7.0" => "570", 
                 "5.7.1" => "571", "6.2.4" => "624", "6.3.1" => "631"}
     vsuffix += "_"+rocmstr[@rocm]


     puts " Do the modules need to be set? (y/n) "
     reply=gets.chomp.strip
     if (reply=="y")

    #    cmds=["module purge",
    #          "module load PrgEnv-amd/8.5.0",
    #          "module load Core",
    #          "module load cmake",
    #          "module load amd/6.1.3",
    #          "module load rocm/6.1.3",
    #          "module load craype-accel-amd-gfx90a",
    #          "module load cray-mpich/8.1.30"]
         if (@env=="cray" and @rocm=="6.3.1")
     #   cmds=["module purge",
     #         "module load Core",
     #         "module load cmake",
     #         "module load PrgEnv-cray",
     #         "module load cpe/24.11",   
     #         "module load rocm/6.3.1",
     #         "module load craype-accel-amd-gfx90a",
     #         "module load cray-hdf5-parallel/1.14.3.3"]

        elsif (@env=="cray" and @rocm=="5.3.0")
        cmds=["module purge",
              "module load Core",
              "module load cmake",
              "module load PrgEnv-cray/8.3.3",
              "module load craype-x86-trento",
              "module load cpe/23.03",   
              "module load rocm/5.3.0",
              "module load craype-accel-amd-gfx90a",
              "module load cray-hdf5-parallel/1.12.2.1"]

        elsif (@env=="amd" and @rocm=="5.3.0")
        cmds=["module purge",
              "module load Core",
              "module load cmake",
              "module load PrgEnv-amd/8.3.3",
              "module load craype-x86-trento",
              "module load amd/5.3.0", 
              "module load rocm/5.3.0",
              "module load cray-mpich/8.1.23",
              "module load cray-libsci/22.12.1.1",
              "module load craype/2.7.31.11",
              "module load craype-accel-amd-gfx90a"]
        elsif (@env=="amd" and @rocm=="5.4.0")
        cmds=["module purge",
              "module load Core",
              "module load cmake",
              "module load PrgEnv-amd/8.5.0",
              "module load craype-x86-trento",
              "module load amd/5.4.0", 
              "module load rocm/5.4.0",
              "module load cray-mpich/8.1.25",
              "module load craype/2.7.31.11",
              "module load cray-libsci/23.12.5",
              "module load craype-accel-amd-gfx90a"]
        elsif (@env=="amd" and @rocm=="5.4.3")
        cmds=["module purge",
              "module load Core",
              "module load cmake",
              "module load PrgEnv-amd/8.5.0",
              "module load craype-x86-trento",
              "module load amd/5.4.3", 
              "module load rocm/5.4.3",
              "module load cray-mpich/8.1.26",
              "module load craype/2.7.31.11",
              "module load cray-libsci/23.12.5",
              "module load craype-accel-amd-gfx90a"]
        elsif (@env=="amd" and @rocm=="5.5.1")
        cmds=["module purge",
              "module load Core",
              "module load cmake",
              "module load PrgEnv-amd/8.5.0",
              "module load craype-x86-trento",
              "module load amd/5.5.1", 
              "module load rocm/5.5.1",
              "module load cray-mpich/8.1.27",
              "module load craype/2.7.31.11",
              "module load cray-libsci/23.12.5",
              "module load craype-accel-amd-gfx90a"]
       # cmds=["module purge",
       #       "module load Core",
       #       "module load cmake",
       #       "module load PrgEnv-amd",
       #       "module load cpe/23.09",
       #       "module load amd/5.5.1", 
       #       "module load rocm/5.5.1",
       #       "module load craype-x86-trento",
       #       "module load craype-accel-amd-gfx90a"]
        elsif (@env=="amd" and @rocm=="5.6.0")
        cmds=["module purge",
              "module load Core",
              "module load cmake",
              "module load PrgEnv-amd/8.5.0",
              "module load craype-x86-trento",
              "module load amd/5.6.0", 
              "module load rocm/5.6.0",
              "module load cray-mpich/8.1.27",
              "module load craype/2.7.31.11",
              "module load cray-libsci/23.12.5",
              "module load craype-accel-amd-gfx90a"]



        elsif (@env=="amd" and @rocm=="5.7.0")
        cmds=["module purge",
              "module load Core",
              "module load cmake",
              "module load PrgEnv-amd/8.5.0",
              "module load craype-x86-trento",
              "module load amd/5.7.0", 
              "module load rocm/5.7.0",
              "module load cray-mpich/8.1.28",
              "module load craype/2.7.31.11",
              "module load cray-libsci/23.12.5",
              "module load craype-accel-amd-gfx90a"]
        elsif (@env=="amd" and @rocm=="5.7.1")
        cmds=["module purge",
              "module load Core",
              "module load cmake",
              "module load PrgEnv-amd/8.5.0",
              "module load craype-x86-trento",
              "module load amd/5.7.1", 
              "module load rocm/5.7.1",
              "module load cray-mpich/8.1.28",
              "module load craype/2.7.31.11",
              "module load cray-libsci/23.12.5",
              "module load craype-accel-amd-gfx90a"]
        elsif (@env=="amd" and @rocm=="6.2.4")
        cmds=["module purge",
              "module load Core",
              "module load cmake",
              "module load PrgEnv-amd/8.6.0",
              "module load craype-x86-trento",
              "module load amd/6.2.4", 
              "module load rocm/6.2.4",
              "module load cray-mpich/8.1.31",
              "module load craype/2.7.31.11",
              "module load cray-libsci/23.12.5",
              "module load craype-accel-amd-gfx90a"]
        elsif (@env=="amd" and @rocm=="6.3.1")
        cmds=["module purge",
              "module load Core",
              "module load cmake",
              "module load PrgEnv-amd/8.6.0",
              "module load craype-x86-trento",
              "module load amd/6.3.1", 
              "module load rocm/6.3.1",
              "module load cray-mpich/8.1.31",
              "module load craype/2.7.31.11",
              "module load cray-libsci/23.12.5",
              "module load craype-accel-amd-gfx90a"]

              #"module load cpe/23.03",   
              #"module load rocm/5.3.0",
              #"module load craype-accel-amd-gfx90a",
              #"module load cray-hdf5-parallel/1.12.2.1"]
         end



     end
  end

      #  load modules if requested:  these are written to a file
      #  and this program exits.  The user must then source
      #  the shell script produced which loads the appropriate modules.
      #  These must be set in the parent process to the process
      #  running this script.

  if (!cmds.empty?)
     env_script=String.new
     env_script+="#!/usr/bin/bash\n"
     env_script+="\n# This file created by build_quda.rb: do not edit\n\n"
     for cmd in cmds
        env_script+=cmd+"\n"
        end
     env_script+="module list\n"
     envfile=File.new("setup_env.sh","w")
     envfile.puts env_script
     envfile.close
     puts "\n\nNecessary module commands written to file setup_env.sh"
     puts "Run \"source setup_env.sh\" to set the environment variables\n\n"
     exit(0)
  end
  
      # default directories

  now=DateTime.now
  nowstr="DATE-"+now.year.to_s+"-"+now.month.to_s+"-"+now.day.to_s

  if (@machine=="cmucluster")
     @qudadir="/home/colin/quda"
     @srcdir=@qudadir+"/source"
     @builddir=@qudadir+"/build/"+nowstr
     @installdir=@qudadir+"/install/"+nowstr
     @qmpdir=""
  elsif (@machine=="perlmutter")
     @qudadir="/global/cfs/cdirs/m2986/stoch_nn/colin/quda"
     @srcdir="/global/cfs/cdirs/m2986/stoch_nn/colin/quda/source"
     @builddir=@qudadir+"/build/"+nowstr
     @installdir="/global/cfs/cdirs/m2986/stoch_nn/colin/quda/install/perlmutter/"+nowstr
     @qmpdir=""
  elsif (@machine=="summit")
     @qudadir="/autofs/nccs-svm1_proj/nph162/colin/quda"
     @srcdir="/autofs/nccs-svm1_proj/nph162/colin/quda/source"
     @builddir=@qudadir+"/build/"+nowstr
     @installdir="/autofs/nccs-svm1_proj/nph162/colin/quda/install/"+nowstr
     @qmpdir="/autofs/nccs-svm1_proj/nph162/colin/usqcd_2023/install/DATE-2024-3-21/qmp_omp/lib/cmake/QMP"
  elsif (@machine=="frontier")
     @qudadir="/autofs/nccs-svm1_proj/nph165/colin/quda"
     @srcdir="/autofs/nccs-svm1_proj/nph165/colin/quda/source"
     @builddir=@qudadir+"/build/"+nowstr
     @installdir="/autofs/nccs-svm1_proj/nph165/colin/quda/install/"+nowstr
     @qmpdir="/autofs/nccs-svm1_proj/nph165/colin/usqcd_2023/install/DATE-2024-3-21/qmp_omp/lib/cmake/QMP"
  elsif (@machine=="bridges")
     @qudadir="/jet/home/mornings/quda"
     @srcdir=@qudadir+"/source"
     @builddir=@qudadir+"/build/"+nowstr
     @installdir=@qudadir+"/install/"+nowstr
     @qmpdir=""
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
  commchoices.push("qmp-parallel")
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
  commstr = { "mpi-parallel" => "mpi", "none-serial" => "ser", "qmp-parallel" => "qmp" }
  suffix = "_"+commstr[@comm]

  if (@comm=="qmp-parallel")
     puts
     puts " The QMP installation cmake directory for this machine is "
     puts
     puts "        "+@qmpdir
     puts
     puts " Hit <enter> to accept the default, or enter a different path: "
     print " QMP installation cmake directory? "
     reply=gets.chomp.strip
     if (reply!="")
        @qmpdir=reply
        puts " The current QMP installation cmake directory is "
        puts "        "+@qmpdir
     end
  end

  @builddir=@builddir+suffix+vsuffix
  @installdir=@installdir+suffix+vsuffix

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
     @gpuextra=""
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
     @gpuextra=" -DCUDA_cublas_LIBRARY=/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/math_libs/12.2/lib64/libcublas.so"+
               " -DCUDA_cufft_LIBRARY=/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/math_libs/12.2/lib64/libcufft.so"
     @makeopt="-j4"
     @cflags="-O3 -fopenmp -target-accel=nvidia80"
     @cxxflags="-O3 -fopenmp -target-accel=nvidia80"
     @linkerflags=""
          
  elsif (@machine=="bridges")
  
     if (@comm=="mpi-parallel")
        #@cc="mpiicc"
        #@cxx="mpiicpc"
        @cc="mpicc"
        @cxx="mpicxx"
     elsif (@comm=="none-serial")
        @cc="icc"
        @cxx="icpc"
     end
     @gpu="cuda"
     @gpudir="/opt/packages/cuda/v12.6.1"
     @gpubindir=@gpudir+"/bin"
     @gpuincdir=@gpudir+"/include"
     @gpuarch="sm_70"
     @gpuextra=""
     @makeopt="-j8"
     @cflags="-O3 -fopenmp"
     @cxxflags="-O3 -fopenmp"
     @linkerflags=""

  elsif (@machine=="summit")
  
     if (@comm=="mpi-parallel")||(@comm=="qmp-parallel")
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
     @gpuextra=""
     @makeopt="-j8"
     @cflags="-O3 -fopenmp"
     @cxxflags="-O3 -fopenmp"
     @linkerflags=" -lsframe"

  elsif (@machine=="frontier")
  
     if (@comm=="mpi-parallel")||(@comm=="qmp-parallel")
        @cc="cc"
        @cxx="CC"
     elsif (@comm=="none-serial")
        abort("Not currently supported")
     end
     @gpu="hip"
     @gpudir="/opt/rocm-"+@rocm
     @gpubindir=@gpudir+"/bin"
     @gpuincdir=@gpudir+"/include"
     @gpuarch="gfx90a"
     @gpuextra=""
     @makeopt="-j8"

  #   @cflags="-O3 -fopenmp -I${MPICH_DIR}/include -Wall --offload-arch=gfx90a"
  #   @cxxflags="-O3 -fopenmp  -I${MPICH_DIR}/include  -Wall --offload-arch=gfx90a"
  #   @linkerflags=" -L${MPICH_DIR}/lib -lmpi ${CRAY_XPMEM_POST_LINK_OPTS} -lxpmem "
  #   @linkerflags+="  ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a}"
  #   @linkerflags+=" --offload-arch=gfx90a"

     @cflags="-O3 -Wall -fopenmp -I${MPICH_DIR}/include -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 --rocm-path=${ROCM_PATH} --offload-arch=gfx90a "
     @cxxflags="-O3 -Wall -fopenmp -I${MPICH_DIR}/include -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 --rocm-path=${ROCM_PATH} --offload-arch=gfx90a "
     @linkerflags=" --rocm-path=${ROCM_PATH} -L${ROCM_PATH}/lib -lamdhip64 "   #    --offload-arch=gfx90a"

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
      qmp="OFF"
   elsif (@comm=="qmp-parallel")
      mpi="OFF"
      qmp="ON"
   else
      mpi="OFF"
      qmp="OFF"
   end

   if (@gpu=="cuda")
      gpuoptions=" -DQUDA_TARGET_TYPE=\"CUDA\""+
                 " -DQUDA_GPU_ARCH="+@gpuarch+
                 " -DCUDAToolkit_BIN_DIR="+@gpubindir
                 " -DCUDAToolkit_CUPTI_INCLUDE_DIR="+@gpuincdir+@gpuextra
   elsif (@gpu=="hip")
      gpuoptions=" -DQUDA_TARGET_TYPE=\"HIP\""+
                 " -DQUDA_GPU_ARCH="+@gpuarch
          #       " -DCMAKE_HIP_COMPILER="+@gpubindir+"/hipcc"
          #       " -DROCM_PATH=\""+@gpudir+"\""+@gpuextra
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
       " -DQUDA_QMP="+qmp+
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
       " -DCMAKE_SHARED_LINKER_FLAGS=\""+@linkerflags+"\""+
       " -DCMAKE_EXE_LINKER_FLAGS=\""+@linkerflags+"\""

   if (@comm=="qmp-parallel")
      cmd+=" -DQMP_DIR=\""+@qmpdir+"\""
   end
   
   puts cmd

     # do cmake
   system(cmd)
   success=($?.exitstatus==0)
   abort("QUDA cmake failed; aborting") if (!success)
     # do make ( due to the "tee", need PIPESTATUS to get exit status of first cmake command)
   system("cmake --build "+@builddir+" "+@makeopt+" 2>&1 | tee "+@builddir+"/build.log; exit ${PIPESTATUS[0]}")
   success=($?.exitstatus==0)
   abort("QUDA build failed; aborting") if (!success)
     # do make install
   system("cmake --install "+@builddir)
   success=($?.exitstatus==0)
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

