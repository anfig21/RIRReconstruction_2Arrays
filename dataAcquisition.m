function Data = dataAcquisition(Data)
%Data = dataAcquisition(Data) Reads the corresponding data of the given
%loudspeaker.
%   Input:
%       - Data      : data structure.
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% SOURCE
Data.Source.pos = h5read(Data.loudspeaker,'/source/position');

%% REFERENCE LINE
% Line 1
Data.Line1.pos = h5read(Data.loudspeaker,'/dataset_lin2/position');
Data.Line1.h = h5read(Data.loudspeaker,'/dataset_lin2/impulse_response');
Data.Line1.n = h5read(Data.loudspeaker,'/dataset_lin2/noise');

% Line 2
Data.Line2.pos = h5read(Data.loudspeaker,'/dataset_lin1/position');
Data.Line2.h = h5read(Data.loudspeaker,'/dataset_lin1/impulse_response');
Data.Line2.n = h5read(Data.loudspeaker,'/dataset_lin1/noise');

% Line 3
Data.Line3.pos = h5read(Data.loudspeaker,'/dataset_lin3/position');
Data.Line3.h = h5read(Data.loudspeaker,'/dataset_lin3/impulse_response');
Data.Line3.n = h5read(Data.loudspeaker,'/dataset_lin3/noise');

%% SPHERICAL ARRAY
% Middle spherical array
% Data.Sph1.pos = h5read(Data.loudspeaker,'/dataset_sph1/position');
% Data.Sph1.h = h5read(Data.loudspeaker,'/dataset_sph1/impulse_response');
% Data.Sph1.p = h5read(Data.loudspeaker,'/dataset_sph1/pressure');
% Data.Sph1.n = h5read(Data.loudspeaker,'/dataset_sph1/noise');
% Data.Sph1.R0 = mean(Data.Sph1.pos,1);

% Left spherical array
Data.SphL.pos = h5read(Data.loudspeaker,'/dataset_sph2/position');
Data.SphL.h = h5read(Data.loudspeaker,'/dataset_sph2/impulse_response');
Data.SphL.p = h5read(Data.loudspeaker,'/dataset_sph2/pressure');
Data.SphL.n = h5read(Data.loudspeaker,'/dataset_sph2/noise');
Data.SphL.R0 = mean(Data.SphL.pos,1);

% Right spherical array
Data.SphR.pos = h5read(Data.loudspeaker,'/dataset_sph3/position');
Data.SphR.h = h5read(Data.loudspeaker,'/dataset_sph3/impulse_response');
Data.SphR.p = h5read(Data.loudspeaker,'/dataset_sph3/pressure');
Data.SphR.n = h5read(Data.loudspeaker,'/dataset_sph3/noise');
Data.SphR.R0 = mean(Data.SphR.pos,1);

disp('Reading data... OK')

end
