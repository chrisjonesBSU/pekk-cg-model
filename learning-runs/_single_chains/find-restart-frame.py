import signac
import numpy as np
import gsd.hoomd
import os

project = signac.get_project()

for kT, jobs in project.find_jobs(
        filter={"tau_p": 0.1}, doc_filter={"done": True}).groupby("kT_quench"):

    job = list(jobs)[0]
    data = np.genfromtxt(job.fn("sim_traj.txt"), names=True)
    vol = data["mdcomputeThermodynamicQuantitiesvolume"][-2000:]
    job.doc.mean_vol_reduced = np.mean(vol)
    job.doc.mean_box_L_reduced = np.mean(vol)**(1/3)

    box_ratios = []
    frames = []
    with gsd.hoomd.open(job.fn("sim_traj.gsd")) as traj:
        for i in range(20):
            snap = traj[-(i+1)]
            box_ratios.append(snap.configuration.box[0]/job.doc.mean_box_L_reduced)
            frames.append(i+1)

    min_idx = np.where(np.array(box_ratios) == min(box_ratios, key = lambda x: abs(x-1)))[0][0]
    job.doc.restart_frame = frames[min_idx]

    with gsd.hoomd.open(job.fn("sim_traj.gsd")) as traj:
        restart_snap = traj[-job.doc.restart_frame]
    with gsd.hoomd.open(os.path.join(job.ws, "npt-restart.gsd"), "wb") as traj2:
        traj2.append(restart_snap)
