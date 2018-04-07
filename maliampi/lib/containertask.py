import luigi
import sciluigi
import logging
import subprocess

# Setup logging
log = logging.getLogger('sciluigi-interface')


class ContainerInfo():
    """
    A data object to store parameters related to running a specific
    tasks in a container (docker / batch / etc)
    """

    # Meant to store tuples of external to internal mounts, as would be done with -v ext:int
    mountpoints = None
    # num vcpu required
    vcpu = None
    # max memory (mb)
    mem = None
    # Env
    env = None
    # Timeout in seconds
    timeout = None
    # engine
    #   Docker by default. Extensible in the future for batch, slurm-singularity, etc

    def __init__(self,
                 mountpoints=[],
                 vcpu=1,
                 mem=4096,
                 env={},
                 timeout=604800):
        self.mountpoints = mountpoints
        self.vcpu = vcpu
        self.mem = mem
        self.env = env
        self.timeout = timeout

    def get_volumes_dict(self, mode='rw'):
        return {mp[0]: {'bind': mp[1], 'mode': mode} for mp in self.mountpoints}

    def __str__(self):
        """
        Return string of this information
        """
        return(
            "cpu {} mem {}".format(
                self.vcpu,
                self.mem,
            ))


class ContainerInfoParameter(sciluigi.parameter.Parameter):
    '''
    A specialized luigi parameter, taking ContainerInfo objects.
    '''
    def parse(self, x):
        if isinstance(x, ContainerInfo):
            return x
        else:
            raise Exception('Parameter is not instance of ContainerInfo: {}. It is instead {}'.format(x, type(x)))


class ContainerHelpers():
    """
    Mixin with various methods and variables for running commands in containers using (Sci)-Luigi
    """
    # Other class-fields
    # Some guidance to whatever eventually runs the container
    # containerinfo = ContainerInfoParameter(default=None)
    # The ID of the container (docker registry style).
    __container__ = None

    def ex(self, command, mounts={}):
        return self.ex_docker(command, mounts)

    def ex_docker(self, command, mounts={}):
        """
        Run command in the container using docker, with mountpoints
        """
        import docker
        client = docker.from_env()
        if isinstance(command, list):
            command = subprocess.list2cmdline(command)
        try:
            log.info("Attempting to run {} in {}".format(
                command,
                self.__container__
            ))
            stdout = client.containers.run(
                image=self.__container__,
                command=command,
                volumes=mounts,
            )
            log.info(stdout)
            return (0, stdout, "")
        except docker.errors.ContainerError:
            log.error("Non-zero return code from the container")
            return (-1, "", "")
        except docker.errors.ImageNotFound:
            log.error("Could not find container {}".format(
                self.__container__)
                )
            return (-1, "", "")
        except docker.errors.APIError:
            log.error("Docker Server failed")
            return (-1, "", "")
        except:
            log.error("Unknown error occurred")
            return (-1, "", "")


# ================================================================================

class ContainerTask(ContainerHelpers, sciluigi.task.Task):
    '''
    luigi task that includes the ContainerHelpers mixin.
    '''
    pass
