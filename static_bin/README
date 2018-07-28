This folder contains statically linked binaries for the following systems:
 - Linux 64bits (built on Ubuntu 16.04 with g++ 5.4.0)
 - Windows 32bits
 - Windows 64bits
The binaries are respectively suffixed with:
 - unix64
 - win32
 - win64

To be able to use the Windows binaries, first identify if you have a 32 or 64
bits system, then remove the associated suffix and add the '.exe' extension.
For the castor-recon program, the final name must be 'castor-recon.exe'.
On linux systems, simply remove the suffix to get the name 'castor-recon'.

Also note that with the 32bits binaries, it is not possible to use datafiles
bigger than 2GB, and all double precision variables inside the CASToR program
are decreased to simple precision so that the benchmark tests will not pass
the cuts. However, the results still remain valid to some extent.

About the configuration folder 'config'. It is also included here. For the
unix binary, the path to the config directory can be set using the environment
variable CASTOR_CONFIG, or using the '-conf' option that overload the value
of the CASTOR_CONFIG variable. For windows binaries, the path of the config
directory is hard-linked to the current difrectory (i.e. as if the config
directory is supposed to be inside the folder where the CASToR program is
executed). However, it can also be changed using the '-conf' option. See the
documentation for further details.
