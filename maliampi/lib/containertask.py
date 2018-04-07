import luigi
import sciluigi
import logging
import subprocess
import docker
import os
from string import Template

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
    container = None
    # Choices include docker right now. Eventually we can add batch, slurm-singularity, etc
    engine = 'docker'


    def map_paths_to_container(self, paths, container_base_path='/mnt'):
        """
        Accepts a dictionary where the keys are identifiers for various targets
        and the value is the HOST path for that target

        What this does is find a common HOST prefix 
        and remaps to the CONTAINER BASE PATH

        Returns a dict of the paths for the targets as they would be seen 
        if the common prefix is mounted within the container at the container_base_path
        """
        common_prefix = os.path.commonprefix(list(paths.values()))
        container_paths = {
            i: os.path.join(
                container_base_path, 
                os.path.relpath(paths[i], common_prefix))
            for i in paths
        }
        return os.path.abspath(common_prefix), container_paths

    def ex(self, command, input_paths={}, output_paths={}, mounts={}):
        if self.engine == 'docker':
            return self.ex_docker(command, input_paths, output_paths, mounts)
        else:
            raise Exception("Container engine {} is invalid".format(self.engine))

    def ex_docker(self, command, input_paths={}, output_paths={}, mounts={}):
        """
        Run command in the container using docker, with mountpoints
        command is assumed to be in python template substitution format
        """
        client = docker.from_env()
        container_paths = {}

        if len(input_paths) > 0:
            input_host_path_ca, input_container_paths = self.map_paths_to_container(
                input_paths,
                container_base_path='/mnt/inputs'
            )
            container_paths.update(input_container_paths)
            mounts[input_host_path_ca]={'bind': '/mnt/inputs', 'mode': 'ro'}
        
        if len(output_paths) > 0:
            output_host_path_ca, output_container_paths = self.map_paths_to_container(
                output_paths,
                container_base_path='/mnt/outputs'
            )
            container_paths.update(output_container_paths)
            mounts[output_host_path_ca]={'bind': '/mnt/outputs', 'mode': 'rw'}

        command = Template(command).substitute(container_paths)

        try:
            log.info("Attempting to run {} in {}".format(
                command,
                self.container
            ))
            stdout = client.containers.run(
                image=self.container,
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
                self.container)
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
