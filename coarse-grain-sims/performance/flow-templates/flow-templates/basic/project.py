"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
import os


class MyProject(FlowProject):
    pass


class Borah(DefaultSlurmEnvironment):
    hostname_pattern = "borah"
    template = "borah.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="shortgpu",
            help="Specify the partition to submit to."
        )


class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="batch",
            help="Specify the partition to submit to."
        )


class Falcon(DefaultSlurmEnvironment):
    hostname_pattern = "ondemand"
    template = "falcon.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="batch",
            help="Specify the partition to submit to."
        )


# Definition of project-related labels (classification)
@MyProject.label
def sim_done(job):
    return job.doc.sim_done


@MyProject.label
def sample_done(job):
    return job.doc.sample_done


@MyProject.post(sim_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="simulate"
)
def run_simulation(job):
    # Add package imports here
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        # Add your script here


@MyProject.pre(sim_done)
@MyProject.post(sample_done)
@MyProject.operation(
        directives={"ngpu": 0, "executable": "python -u"}, name="sample"
)
def sample(job):
    # Add package imports here
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        # Add your script here


if __name__ == "__main__":
    MyProject().main()
