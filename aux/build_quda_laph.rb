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

  puts
  puts " QUDA_LAPH software installation on machine '"+@machine+"'"
  puts
  cmds=[]
  
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

     puts " Do the modules need to be set? (y/n) "
     reply=gets.chomp.strip
     if (reply=="y")
        cmds=["module purge",
              "module load PrgEnv-gnu/8.5.0",
              "module load gpu",
              "module load cmake/3.24.3",
              "module load cpe-cuda/23.12",
              "module load cudatoolkit/12.2",
              "module load craype-accel-nvidia80",
              "module load craype-x86-milan",
              "module load eigen/3.4.0"]
     end
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

  if (@machine=="frontier")

     puts " Do the modules need to be set? (y/n) "
     reply=gets.chomp.strip
     if (reply=="y")
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
              "module load craype-accel-amd-gfx90a",
              "module load openblas/0.3.26-omp"]
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
     @qudalaphdir="/home/colin/quda_laph"
     @srcdir=@qudalaphdir+"/source"
     @builddir=@qudalaphdir+"/build/"+nowstr
     @installdir=@qudalaphdir+"/install/"+nowstr
     @qudainstalldir="/home/colin/quda/install/DATE-2025-1-9"
  elsif (@machine=="perlmutter")
     @qudalaphdir="/global/cfs/cdirs/m2986/stoch_nn/colin/quda_laph"
     @srcdir=@qudalaphdir+"/source"
     @builddir=@qudalaphdir+"/build/"+nowstr
     @installdir=@qudalaphdir+"/install/"+nowstr
     @qudainstalldir="/global/cfs/cdirs/m2986/stoch_nn/colin/quda/install/perlmutter/DATE-2024-11-7"
  elsif (@machine=="summit")
     @qudalaphdir="/autofs/nccs-svm1_proj/nph162/colin/quda_laph"
     @srcdir=@qudalaphdir+"/source"
     @builddir=@qudalaphdir+"/build/"+nowstr
     @installdir=@qudalaphdir+"/install/"+nowstr
     @qudainstalldir="/autofs/nccs-svm1_proj/nph162/colin/quda/install/DATE-2024-11-5"
  elsif (@machine=="frontier")
     @qudalaphdir="/autofs/nccs-svm1_proj/nph165/colin/quda_laph"
     @srcdir=@qudalaphdir+"/source"
     @builddir=@qudalaphdir+"/build/"+nowstr
     @installdir=@qudalaphdir+"/install/"+nowstr
     @qudainstalldir="/autofs/nccs-svm1_proj/nph165/colin/quda/install/DATE-2025-1-24"
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
     @hostcblas="openblas"
     @cblaslib="/opt/cray/pe/libsci/23.12.5/GNU/12.3/x86_64/lib/libsci_gnu_mpi.so"
     @cblasinc="/opt/cray/pe/libsci/23.12.5/GNU/12.3/x86_64/include"
          
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
     #@hostcblas="gsl"
     #@cblaslib="/sw/summit/spack-envs/summit-plus/opt/gcc-8.5.0/gsl-2.7.1-643fyt4m7sv5d4qjsstofh2l5rlqly2j/lib/libgslcblas.so"
     #@cblasinc="/sw/summit/spack-envs/summit-plus/opt/gcc-8.5.0/gsl-2.7.1-643fyt4m7sv5d4qjsstofh2l5rlqly2j/include"
     @hostcblas="openblas"
     @cblaslib="/sw/summit/spack-envs/summit-plus/opt/gcc-12.1.0/openblas-0.3.24-hwvdhvhbfg4odxpj7mhv6ij3jtitcmyt/lib/libopenblas.so"
     @cblasinc="/sw/summit/spack-envs/summit-plus/opt/gcc-12.1.0/openblas-0.3.24-hwvdhvhbfg4odxpj7mhv6ij3jtitcmyt/include"

  elsif (@machine=="frontier")
  
     if (@comm=="mpi-parallel")||(@comm=="qmp-parallel")
        @cc="cc"
        @cxx="CC"
     elsif (@comm=="none-serial")
        abort("Not currently supported")
     end
     @gpu="hip"
     @gpudir="/opt/rocm-6.3.1"
     #@gpubindir=@gpudir+"/bin"
     #@gpuincdir=@gpudir+"/include"
     #@gpuarch="gfx90a"
     @makeopt="-j32"
     @cflags=" -D__HIP_PLATFORM_AMD__ -O2  -fopenmp"
     @cxxflags=" -D__HIP_PLATFORM_AMD__ -O2 -fopenmp"
     @linkerflags=""
     @hostcblas="openblas"
     @cblaslib="/sw/frontier/spack-envs/core-24.07/opt/gcc-7.5.0/openblas-0.3.26-2uukuglbn3pfur3tzlfe5dxmj6fjauoi/lib/libopenblas.so"
     @cblasinc="/sw/frontier/spack-envs/core-24.07/opt/gcc-7.5.0/openblas-0.3.26-2uukuglbn3pfur3tzlfe5dxmj6fjauoi/include"
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

