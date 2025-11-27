from argparse import ArgumentParser
from pathlib import Path
from nolanlab_ephys.si_protocols import protocols
from nolanlab_ephys.utils import get_recording_folders, chronologize_paths
from nolanlab_ephys.probe_info import make_probe_plot
from nolanlab_ephys.si_protocols import generic_postprocessing, compute_automated_curation

import spikeinterface.full as si

def main():

    parser = ArgumentParser()

    parser.add_argument('mouse')
    parser.add_argument('day')
    parser.add_argument('sessions')
    parser.add_argument('protocols')
    parser.add_argument('--data_folder', default=None)
    parser.add_argument('--deriv_folder', default=None)

    mouse = int(parser.parse_args().mouse)
    day = int(parser.parse_args().day)

    sessions_string = parser.parse_args().sessions
    sessions = sessions_string.split(',')

    protocols = parser.parse_args().protocols
    protocols_list = protocols.split(',')

    data_folder = parser.parse_args().data_folder
    if data_folder is None:
        data_folder = "/home/nolanlab/Work/Harry_Project/data/"
    data_folder = Path(data_folder)

    deriv_folder = parser.parse_args().deriv_folder
    if deriv_folder is None:
        deriv_folder = "/home/nolanlab/Work/Harry_Project/derivatives/"
    deriv_folder = Path(deriv_folder)

    mouseday_deriv_folder = deriv_folder / f"M{mouse:02d}/D{day:02d}"
    mouseday_deriv_folder.mkdir(parents=True, exist_ok=True)

    si.set_global_job_kwargs(n_jobs=8)

    for protocol in protocols_list:
        _ = do_sorting_pipeline(mouse, day, sessions, data_folder, deriv_folder, protocol)

def do_sorting_pipeline(mouse, day, sessions, data_folder, deriv_folder, protocol):

    protocol_info = protocols[protocol]

    recording_paths = chronologize_paths(get_recording_folders(data_folder=data_folder, mouse =mouse, day=day))

    make_probe_plot(recording_paths[0], save_path=deriv_folder / f"M{mouse:02d}/D{day:02d}/M{mouse:02d}_D{day:02d}_probe_layout.png")

    recordings = [si.read_openephys(recording_path) for recording_path in recording_paths]

    concatenated_recording = si.concatenate_recordings(recordings)

    pp_recording = si.apply_preprocessing_pipeline(concatenated_recording, protocol_info['preprocessing'])

    sorting = si.run_sorter(recording=pp_recording, **protocol_info['sorting'], remove_existing_folder=True, verbose=True, folder=f"M{mouse:02d}_D{day:02d}_{protocol}_{'-'.join(sessions)}_output")

    for recording, session in zip(recordings, sessions):

        # we do all our syncing assuming that t=0 is at the start of the ephys data
        recording._recording_segments[0].t_start = 0

        analyzer = si.create_sorting_analyzer(
            recording=si.apply_preprocessing_pipeline(recording, protocol_info['preprocessing_for_analyzer']), 
            sorting=sorting, 
            folder = deriv_folder / f"M{mouse:02d}/D{day:02d}/{session}/{protocol}/sub-{mouse:02d}_day-{day:02d}_ses-{session}_srt-{protocol}_analyzer",
            format = "zarr",
            peak_sign = "both",
            radius_um = 70,
        )

        analyzer.compute(generic_postprocessing)

    unitrefine_model_path = "/exports/eddie/scratch/chalcrow/wolf/code/models/curationA_model"
    curation_output_path = deriv_folder / f"M{mouse:02d}/D{day:02d}/{session}/{protocol}/sub-{mouse:02d}_day-{day:02d}_ses-{session}_srt-{protocol}_curationA.json"
    compute_automated_curation(analyzer, unitrefine_model_path, curation_output_path=curation_output_path)

    return analyzer

if __name__ == "__main__":
    main()
