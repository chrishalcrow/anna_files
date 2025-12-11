from neuroconv.tools.spikeinterface import add_sorting_to_nwbfile
import spikeinterface.full as si
import pandas as pd
import pynwb
import numpy as np

metadata_path = derivatives_path / "metadata" / "generic_metadata.yml"
nwbfile = pynwb.NWBFile(**metadata["NWBFile"])

#########
# UNITS #
#########

# spikeinterface stuff

sorting_analyzer = si.load_sorting_analyzer(analyzer_path)
sorting = sorting_analyzer.sorting

add_sorting_to_nwbfile(sorting, nwbfile)

# we store the...

# quality metrics
quality_metrics = sorting_analyzer.get_extension('quality_metrics').get_data()

# template metrics
template_metrics = sorting_analyzer.get_extension('template_metrics').get_data()

# unit locations
probe_locations = sorting_analyzer.get_extension("unit_locations").get_data()

probe_locations_df = pd.DataFrame()
for coord, data in zip(["x", "y", "z"], probe_locations.T, strict=True):
    probe_locations_df[f"coord_est_{coord}"] = data

# and which channel the extremal template (template with max amp) is on
probe_locations_df["extremum_channel"] = si.get_template_extremum_channel(
    sorting_analyzer
)

all_unit_metadata = pd.concat([probe_locations_df, quality_metrics, template_metrics], axis=1)

# here, we add all that to the nwbfile
for column_name, unit_data in all_unit_metadata.items():
    nwbfile.units.add_column(name=column_name, data=unit_data.to_numpy(), description=f"{column_name}, computed using spikeinterface.")

# Then we also have all the locations in CCFs coords, which we've already saved in a csv
all_cluster_annotations_path = derivatives_path / "labels" / "anatomy" / "cluster_annotations.csv"
all_cluster_annotations = pd.read_csv(all_cluster_annotations_path)
cluster_annotations = all_cluster_annotations.query(f'mouse == {mouse} & day == {day}')

unit_coord_data = cluster_annotations[['coord_CCFs_z', 'coord_CCFs_y', 'coord_CCFs_x']].values
unit_region_data = np.array(cluster_annotations['brain_region'].values, dtype=str)

nwbfile.units.add_column(data=unit_coord_data[:,0], name="coord_CCFs_z", description="z coordinate in Common Coordinate Framework")
nwbfile.units.add_column(data=unit_coord_data[:,1], name="coord_CCFs_y", description="y coordinate in Common Coordinate Framework")
nwbfile.units.add_column(data=unit_coord_data[:,2], name="coord_CCFs_x", description="x coordinate in Common Coordinate Framework")
nwbfile.units.add_column(data=unit_region_data, name="brain_region", description="Brain region according to Allen Atlas")
