# Setup to run 64-bit MATLAB executables on LSC grid computers
#   MATLAB version R2008a
#
export MATLAB_ROOT=/ldcg/matlab_r2013a
export ARCH=glnxa64

# Add MCR location in MATLAB distrib to LD_LIBRARY_PATH
#   We check that LD_LIBRARY_PATH is defined and that we
#       haven't already added MCR to this
#    ** we add sys/opengl/lib/glnxa64 to cover dependencies missed
#           in previous compiler revs
#    ** we have added /bin/glnxa64 due to some known compiler problems if
#               tracking down all dependencies
#
if [ ! "$LD_LIBRARY_PATH" ] 
  then
    export LD_LIBRARY_PATH=${MATLAB_ROOT}/sys/opengl/lib/glnxa64
    export LD_LIBRARY_PATH=${MATLAB_ROOT}/sys/java/jre/glnxa64/jre1.6.0/lib/amd64:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=${MATLAB_ROOT}/sys/java/jre/glnxa64/jre1.6.0/lib/amd64/server:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=${MATLAB_ROOT}/sys/java/jre/glnxa64/jre1.6.0/lib/amd64/native_threads:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=${MATLAB_ROOT}/extern/lib/glnxa64:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=${MATLAB_ROOT}/bin/glnxa64:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=${MATLAB_ROOT}/sys/os/glnxa64:${LD_LIBRARY_PATH}
  else
    echo $LD_LIBRARY_PATH | grep ${MATLAB_ROOT}/runtime/glnxa64 > /dev/null
    if [ "$?" -ne 0 ]  
      then
        export LD_LIBRARY_PATH=${MATLAB_ROOT}/sys/opengl/lib/glnxa64:${LD_LIBRARY_PATH}
        export LD_LIBRARY_PATH=${MATLAB_ROOT}/sys/java/jre/glnxa64/jre1.6.0/lib/amd64:${LD_LIBRARY_PATH}
        export LD_LIBRARY_PATH=${MATLAB_ROOT}/sys/java/jre/glnxa64/jre1.6.0/lib/amd64/server:${LD_LIBRARY_PATH}
        export LD_LIBRARY_PATH=${MATLAB_ROOT}/sys/java/jre/glnxa64/jre1.6.0/lib/amd64/native_threads:${LD_LIBRARY_PATH}
        export LD_LIBRARY_PATH=${MATLAB_ROOT}/extern/lib/glnxa64:${LD_LIBRARY_PATH}
        export LD_LIBRARY_PATH=${MATLAB_ROOT}/bin/glnxa64:${LD_LIBRARY_PATH}
        export LD_LIBRARY_PATH=${MATLAB_ROOT}/sys/os/glnxa64:${LD_LIBRARY_PATH}
    fi
fi
#export XAPPLRESDIR=${MATLAB_ROOT}/X11/app-defaults
export MCR_CACHE_ROOT=/usr1/${USER}
# export MCR_CACHE_ROOT=/usr/local/MCR/matlab_r2013a/v81

#MCRROOT=/usr/local/MCR/matlab_r2013a/v81;

MCRROOT=/ldcg/matlab_r2013a;

#MCRROOT=/ldcg/matlab_r2013a/bin/glnxa64;

export XAPPLRESDIR=${MCRROOT}/X11/app-defaults;

export LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64;
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64 ;
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ;
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;

# clear any setting of MATLAB, which is set when MATLAB is started
export -n MATLAB
