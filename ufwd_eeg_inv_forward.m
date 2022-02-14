function D = ufwd_eeg_inv_forward(varargin)
% Compute M/EEG leadfield with forward models currently unsupported with
% FieldTrip
% FORMAT D = spm_eeg_inv_forward(D,val)
%
% D                - input struct
% (optional) fields of S:
% D                - filename of EEG/MEG mat-file
%
% Output:
% D                - EEG/MEG struct with filenames of Gain matrices)
%__________________________________________________________________________
% Copyright (C) 2008-2021 Wellcome Centre for Neuroimaging

% George O'Neill, Jeremie Mattout & Christophe Phillips
% $Id$


%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename);

%-Initialisation
%--------------------------------------------------------------------------
[D, val] = spm_eeg_inv_check(varargin{:});

if numel(D.inv{val}.datareg) ~= numel(D.inv{val}.forward)
    error('Separate coregistration is required for every modality.');
end

Fgraph = spm_figure('FindWin','Graphics');
spm_figure('Clear',Fgraph);
spm('Pointer', 'Watch');
if isempty(Fgraph) || spm('CmdLine'), graph = 'no'; else, graph = 'yes'; end

for i = 1:numel(D.inv{val}.forward)
    M    = D.inv{val}.datareg(i).fromMNI*D.inv{val}.mesh.Affine;
    
    M    = diag([1e-3 1e-3 1e-3 1])*M; % convert to m
    
    mesh = spm_eeg_inv_transform_mesh(M, D.inv{val}.mesh);
    
    mesh_correction = [];
    
    sens = D.inv{val}.datareg(i).sensors;
    
    siunits = true;
    if isequal(D.inv{val}.datareg(i).modality, 'MEG')
        if  isequal(unique(sens.chanunit(strmatch('MEG', sens.chantype))), {'snr'})
            sens = ft_datatype_sens(sens, 'distance', 'm');
            siunits = false;
        else
            sens = ft_datatype_sens(sens, 'amplitude', 'T', 'distance', 'm');
        end
    else
        if  isequal(unique(sens.chanunit(strmatch('EEG', sens.chantype))), {'snr'})
            sens = ft_datatype_sens(sens, 'distance', 'm');
            siunits = false;
        else
            sens = ft_datatype_sens(sens, 'amplitude', 'V', 'distance', 'm');
        end
    end
    
    switch D.inv{val}.forward(i).voltype
        case 'uBEM'
            
            % if using a template MRI, save the BEM file to the directory
            % of the MEEG mat/dat files, otherwise save next to MRI.
            if mesh.template
                [~,name,~,~] = spm_fileparts(mesh.sMRI);
                volfile = fullfile(D.path,[name '_uBEM.mat']);
            else
                volfile = spm_file(mesh.sMRI, 'suffix','_uBEM', 'ext','mat');
            end
            vol = [];
            
            if exist(volfile, 'file')
                vol = ft_read_headmodel(volfile);
                if ~isfield(vol, 'unit') || ~isequal(vol.unit, 'm')
                    vol = [];
                end
            end
            
            if isempty(vol)
                vol.cond   = [0.3300 0.0041 0.3300];
                % brain
                vol.bnd(1) = export(gifti(mesh.tess_iskull), 'ft');
                % skull
                vol.bnd(2) = export(gifti(mesh.tess_oskull), 'ft');
                % skin
                vol.bnd(3) = export(gifti(mesh.tess_scalp),  'ft');
                
                cfg                         = [];
                cfg.method                 = 'ubem';
                cfg.siunits                = 'yes';
                cfg.showcallinfo           = 'no';
                vol                        = ufwd_prepare_headmodel(cfg, vol);
                save(volfile, 'vol', spm_get_defaults('mat.format'));
                
            end
            
            cfg = [];
            cfg.headmodel = vol;
            cfg.grid.pos = mesh.tess_ctx.vert;
            cfg.moveinward = 6e-3; %move to empirically determined BEM safe zone
            gridcorrect = ufwd_prepare_sourcemodel(cfg);
            
            mesh_correction    = rmfield(cfg, {'headmodel', 'grid'});
            
            mesh.tess_ctx.vert = gridcorrect.pos;
            
            vol = volfile;
            
            modality            = D.inv{val}.datareg(i).modality;
            
        otherwise
            error('Unsupported volume model type.');
    end
    
    D.inv{val}.forward(i).vol             = vol;
    D.inv{val}.forward(i).mesh            = mesh.tess_ctx;
    D.inv{val}.forward(i).mesh_correction = mesh_correction;
    D.inv{val}.forward(i).modality        = modality;
    D.inv{val}.forward(i).siunits         = siunits;
    
    D.inv{val}.forward(i).sensors  = sens;
    
    D.inv{val}.forward(i).toMNI    = D.inv{val}.datareg(i).toMNI*diag([1e3 1e3 1e3 1]);
    D.inv{val}.forward(i).fromMNI  = diag([1e-3 1e-3 1e-3 1])*D.inv{val}.datareg(i).fromMNI;
    
    spm_figure('Clear',Fgraph);
end

% This is to force recomputing the lead fields
try, D.inv{val} = rmfield(D.inv{val}, 'gainmat'); end

fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('Pointer', 'Arrow');
