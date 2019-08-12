from . import extract_sim
from .. import utils

import os
import time
import logging

LOG = logging.getLogger(__name__)

# JWST imports
from mirisim import MiriSimulation
from mirisim import config_parser

__all__ = ["run"]


def run(observations, scene, simulator=None, dryrun=False):
    """
    Given the list of observations retrieved from the APT file, run the corresponding MIRISim simulations

    For each separate simulation, if the folder she's supposed to end into already exist, the simulation is skipped.
    That means if you don't delete all the output folders from one run to the others, you will not run anything. This
    is intended when something fail, to avoid re-running all the simulations while only the last 3 failed.

    :param observations: list of observations dictionnaries
    :param SceneConfig scene: scene needed for the observations. Each scene must be a SceneConfig object
                  If given only one scene, it will be used for all observations.
                  Else, you need to provide a dictionnary with key being the name of each target in the APT file,
                  and the value being the corresponding SceneConfig object. For instance, in the APT you have
                  "1 NGC-188-FTS107-1", the key will have to be "1 NGC-188-FTS107-1" (including the number at the start)

    :param SimulatorConfig simulator: [optional] if given, replace the default simulator.ini input file
    :param bool dryrun: If true, do not run the simulations but create the corresponding folders/files
    :return:
    """

    simulations = extract_sim.get_simulations(observations)
    nb_simulations = len(simulations)

    # We define the Simulation parameters for MIRISIM
    if simulator is None:
        simulator_config = config_parser.SimulatorConfig.from_default()
    else:
        simulator_config = simulator

    nb_skipped = 0
    start_time = time.time()
    for (sim_id, sim) in enumerate(simulations):
        sim_start = time.time()

        ellapsed = sim_start - start_time
        eta = ellapsed * (nb_simulations - nb_skipped) / (sim_id + 1 - nb_skipped) - ellapsed
        eta_str = utils.strtime(eta)
        msg = "Simulation {}/{} ({}) ; ETA: {}      ".format(sim_id+1, nb_simulations, sim["name"], eta_str)

        print(msg, end="\r", flush=True)

        sim_config = config_parser.SimConfig.makeSim(**sim)

        sim_dir = sim["name"]

        # If the folder exists, we skip this simulation
        if os.path.isdir(sim_dir):
            LOG.info("{} already exists, skipping.".format(sim_dir))
            nb_skipped += 1
            continue

        # If one scene is given, used for all simulations. Else, each simulation must have its scene_file
        if isinstance(scene, config_parser.SceneConfig):
            scene_config = scene
        elif isinstance(scene, dict):
            obs = observations[sim_id]
            target_name = obs["target_id"]
            try:
                scene_config = scene[target_name]
            except KeyError:
                print("\n")  # To leave alone the "progressbar" line.
                LOG.exception("Missing target in input scene dict: {}".format(target_name))
                raise
        else:
            print("\n")  # To leave alone the "progressbar" line.
            raise ValueError("scene must be a SceneConfig, or a dict of {{target_id: SceneConfig}}")

        # We prepare the simulation.
        # NOTE: We want to create it during a dryrun too, to make sure everything would work in a real script
        # loglevel : "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"
        ms = MiriSimulation(sim_config=sim_config, scene_config=scene_config, simulator_config=simulator_config,
                            loglevel='ERROR')

        if dryrun:
            os.mkdir(sim_dir)
            sim_config.write(filename=os.path.join(sim_dir, 'simulation.ini'))
            scene_config.write(filename=os.path.join(sim_dir, 'scene.ini'))
            simulator_config.write(filename=os.path.join(sim_dir, 'simulator.ini'))
        else:
            ms.run()

            # Renaming the folder
            os.rename(ms.path_out, sim_dir)

        sim_time = utils.strtime(time.time()-sim_start)
        LOG.info("{} completed in {}".format(sim["name"], sim_time))

    completion_time = utils.strtime(time.time()-start_time)
    print("Ran {} simulations, {} skipped, ellapsed time: {}".format(nb_simulations, nb_skipped, completion_time))