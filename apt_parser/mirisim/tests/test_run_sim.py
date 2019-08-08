import pytest

sim = pytest.importorskip("apt_parser.mirisim")
from apt_parser.mirisim import run
import mirisim.config_parser as cp
import os
import shutil

def clean_folders(observations):
    """
    After a test, needs to clean the folders created by the test

    :param observations:
    :return:
    """

    for obs in observations:
        folder = obs["name"]

        if os.path.isdir(folder):
            print("rm {}".format(folder))
            shutil.rmtree(folder)



def test_mirisim_run():
    """"
    Test apt_parser.mirisim.run
    """

    observations = [{'template':"MiriImaging", "filter":"F560W", "obs_id":1, "subarray":"FULL",
               "exposures":1, "ima_frames":1, "ima_integrations":1, "NDither":1,
            "readDetect":"FAST", "name":"obs1_IMA_F560W"}]

    default_scene = cp.SceneConfig.from_default()
    default_scene["sky"]["name"] = "test1"

    scene = {"my_target1":default_scene}

    # Clean any folder from previous failed test
    clean_folders(observations)

    with pytest.raises(ValueError):
        run(observations, scene="test.ini", simulator=None, dryrun=True)

    clean_folders(observations)

    with pytest.raises(KeyError):
        observations[0]["target_id"] = "my_target2"

        run(observations, scene=scene, simulator=None, dryrun=True)

    clean_folders(observations)

    observations.append({'template':"MiriImaging", "filter":"F560W", "obs_id":2, "subarray":"FULL",
               "exposures":1, "ima_frames":1, "ima_integrations":1, "NDither":1,
            "readDetect":"FAST", "name":"obs2_IMA_F560W", "target_id":"my_target1"})

    ref_scene_2 = cp.SceneConfig.from_default()
    ref_scene_2["sky"]["name"] = "test2"
    scene["my_target2"] = ref_scene_2

    # Run a different scene per observation
    run(observations, scene=scene, simulator=None, dryrun=True)
    for obs in observations:
        folder = obs["name"]
        ref_scene = scene[obs["target_id"]]
        test_scene = cp.SceneConfig(os.path.join(folder, "scene.ini"))
        assert ref_scene["sky"]["name"] == test_scene["sky"]["name"], "Wrong SceneConfig used in apt_parser.mirisim.run for {}".format(folder)

    clean_folders(observations)

    # Run the same scene for all observations
    run(observations, scene=default_scene, simulator=None, dryrun=True)
    for obs in observations:
        folder = obs["name"]
        test_scene = cp.SceneConfig(os.path.join(folder, "scene.ini"))
        assert default_scene["sky"]["name"] == test_scene["sky"]["name"], "Wrong SceneConfig used in apt_parser.mirisim.run for {}".format(folder)
    clean_folders(observations)