from pathlib import Path

import spikeinterface.full as si
from spikeinterface.curation.curation_tools import resolve_merging_graph
from spikeinterface.curation import validate_curation_dict

protocols = {

    'kilosort4A': {
        'preprocessing': {
            'detect_and_remove_bad_channels': {'seed': 1205},
            'phase_shift': {},
        },
        'sorting': {
            'sorter_name': 'kilosort4',
            'do_correction': False,
            'use_binary_file': False,
        },
        'preprocessing_for_analyzer': {
            'detect_and_remove_bad_channels': {'seed': 1205},
            'phase_shift': {},
            'common_reference': {},
            'bandpass_filter': {},
        },
    },

    'mountainsort5A': {
        'preprocessing': {
            # 'bandpass_filter': {},
            # 'whiten': {'dtype': 'float'},
        },
        'sorting': {
            'sorter_name': 'mountainsort5',
            'scheme': '2',
        },
        'preprocessing_for_analyzer': {
            'common_reference': {},
            'bandpass_filter': {},
        },
    },

    'mountainsort5B': {
        'preprocessing': {
            # 'bandpass_filter': {},
            # 'whiten': {'dtype': 'float'},
        },
        'sorting': {
            'sorter_name': 'mountainsort5',
            'scheme': '3',
        },
        'preprocessing_for_analyzer': {
            'common_reference': {},
            'bandpass_filter': {},
        },
    },

    'mountainsort4A': {
        'preprocessing': {
        },
        'sorting': {
            'sorter_name': 'mountainsort4',
        },
        'preprocessing_for_analyzer': {
            'common_reference': {},
            'bandpass_filter': {},
        },
    },

    'herdingspikesA': {
        'preprocessing': {
            'bandpass_filter': {},
            'common_reference': {},
     #       'whiten': {'dtype': 'float'},
        },
        'sorting': {
            'sorter_name': 'herdingspikes',
        },
        'preprocessing_for_analyzer': {
            'bandpass_filter': {},
            'common_reference': {},
        },
    },

    'spykingcircus2A': {
        'preprocessing': {
        },
        'sorting': {
            'sorter_name': 'spykingcircus2',
            'apply_motion_correction': False,
        },
        'preprocessing_for_analyzer': {
            'bandpass_filter': {},
            'common_reference': {},
        },
    },

    'tridesclous2A': {
        'preprocessing': {
        },
        'sorting': {
            'sorter_name': 'tridesclous2',
        },
        'preprocessing_for_analyzer': {
            'bandpass_filter': {},
            'common_reference': {},
        },
    },


}

generic_postprocessing = {
    'unit_locations': {},
    'random_spikes': {},
    'noise_levels': {},
    'waveforms': {},
    'templates': {},
    'spike_amplitudes': {'peak_sign': 'both'},
    'isi_histograms': {},
    'spike_locations': {'spike_retriver_kwargs': {'peak_sign': 'both'}},
    'correlograms': {},
    'template_similarity': {'method': 'l2'},
    'quality_metrics': {'metric_names': ['num_spikes', 'firing_rate', 'presence_ratio', 'snr', 'isi_violation', 'rp_violation', 'sliding_rp_violation', 'amplitude_cutoff', 'amplitude_median', 'amplitude_cv', 'synchrony', 'firing_range', 'drift', 'sd_ratio'], 'metric_params': {'snr': {'peak_sign': 'both'}, 'amplitude_cutoff': {'peak_sign': 'both'}, 'amplitude_median': {'peak_sign': 'both'}}},
    'template_metrics': {'include_multi_channel_metrics': True, 'peak_sign': 'both'},
}


def compute_automated_curation(analyzer, model_path, curation_output_path):

    unitrefine_labels = si.auto_label_units(analyzer, model_folder=model_path, trust_model=True)

    noise_units = list(unitrefine_labels[unitrefine_labels['prediction'] == 'noise'].index)

    label_definitions = {
        "quality": dict(name="quality", label_options=["good", "noise"], exclusive=True),
    }

    manual_labels = [
        {"unit_id": unit_id, "labels": {"quality": [unitrefine_labels.loc[unit_index].values[0]]}}
        for unit_index, unit_id in enumerate(analyzer.unit_ids)
    ]

    potential_merges = si.compute_merge_unit_groups(
        analyzer.select_units(list(unitrefine_labels.query('prediction == "good"').index)),
        preset='x_contaminations'
    )
    more_potential_merges = si.compute_merge_unit_groups(
        analyzer.select_units(list(unitrefine_labels.query('prediction == "good"').index)),
        preset="temporal_splits"
    )

    final_merges = resolve_merging_graph(analyzer.sorting, potential_merges + more_potential_merges)

    curation_dict = dict(
        format_version="2",
        unit_ids=analyzer.unit_ids,
        label_definitions=label_definitions, 
        merges=[
            dict(unit_ids=p)
            for p in final_merges
        ],
        removed=noise_units,
        manual_labels=manual_labels
    )

    validate_curation_dict(curation_dict)

    # we can use the CurationModel for more secure serialization
    from spikeinterface.curation.curation_model import CurationModel

    curation = CurationModel(**curation_dict)
    curation_file = Path(curation_output_path)
    _ = curation_file.write_text(curation.model_dump_json(indent=4))
