clear; clc;
% Start timer
tic;



% Example Usage: simple common-source

% Create SCA instance
sca = SCA_new();
sca = sca.set_debug(sca, true);

% Build a simple common-source amplifier
sca = sca.add_mos(sca, 'M1', {'vout', 'vin', '0', '0'}, struct('gm',sym('gm_in'), 'ro',sym('ro_in'), 'cgs',sym('Cgs_in'),'cgd',sym('0')));
% sca = sca.add_mos(sca, 'M2', {'vout', 'v_bias', '0', 'vdd'}, struct('gm',sym('gm_tail'), 'ro',sym('ro_tail'), 'cgs',sym('Cgs_tail'),'cgd',sym('Cgd_tail')));
% sca = add_cap(sca, "Cgs", {'vin', '0'}, sym('Cgs'));
% sca = add_cap(sca, "Cgd", {'vin', 'vout'}, sym('Cgd'));

sca = sca.add_res(sca, 'Rsig', {'vsig', 'vin'}, 'Rsig');
sca = sca.add_res(sca, 'RL', {'vout', 'vdd'}, 'RL');
sca = sca.add_cap(sca, 'CL', {'vout', '0'}, 'CL');

sca = sca.add_cap(sca, 'CC', {'vout', 'vin'}, 'CC');

% Define node role
sca = sca.define_node_role(sca, 'vdd','supply');

% Define I/O
sca = sca.define_io(sca, {'vsig'}, {'vout'}, 'vin');

% Get transfer function
[H, sca] = sca.get_tf(sca);


% Get DC gain
Hdc = sca.get_dc_gain(sca, H);
Hdc_parallel = replace_parallel(Hdc);
pretty(Hdc_parallel);
% % Get poles and zeros
% [poles, zeros] = sca.get_poles_zeros(sca, H);
% 
% Rin_sig = sca.res_looking_from_node(sca, 'vin');
% 
% Cin_Vin = sca.cap_looking_from_node(sca, 'vin');
% 
% Cin_Vout = sca.cap_looking_from_node(sca, 'vout');

% Show netlist
sca.show_netlist(sca);
sca.validate_netlist(sca);



%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Create SCA instance
% sca = SCA_new();
% sca = sca.set_debug(sca, true);
%
% % Parameters
% syms Rsig RL ro gm gmb Cgs Cgd Cgb s
%
% % Signal source resistance
% sca = sca.add_res(sca, 'Rsig', {'vsig', 'G'}, Rsig);
%
% sca = sca.add_res(sca, 'RL', {'D', 'VDD'}, RL);
% sca = sca.add_res(sca, 'ro', {'D', 'S'}, ro);
%
%
% % Gate capacitances
% sca = sca.add_cap(sca, 'Cgs', {'G', 'S'}, Cgs);
% sca = sca.add_cap(sca, 'Cgd', {'G', 'D'}, Cgd);
% % sca = sca.add_cap(sca, 'Cgb', {'G', 'B'}, Cgb);
%
% % --- gm*vgs : VCCS from D->S controlled by V(G,S) ---
% add_vccs(sca, 'gm_vgs', {'D','S'}, {'G','S'}, gm);
%
% % --- gmb*vbs : VCCS from D->S controlled by V(B,S) ---
% add_vccs(sca, 'gmb_vbs', {'D','S'}, {'D','S'}, gmb);
%
% % Ground source, drain, bulk for simplicity
% % sca = sca.add_res(sca, 'R_s', {'S','0'}, 0);  % AC ground
% % sca = sca.add_res(sca, 'R_d', {'D','0'}, 0);
% % sca = sca.add_res(sca, 'R_b', {'B','0'}, 0);
%
% sca = sca.define_node_role(sca,'VDD','supply');
% sca = sca.define_node_role(sca,'S','ground');
% sca = sca.define_node_role(sca,'B','ground');
%
% % Define I/O
% sca = sca.define_io(sca, {'vsig'}, {'D'}, 'vin');
%
% % Get transfer function
% H = sca.get_tf(sca, true);
%
%
% % Get DC gain
% Hdc = sca.get_dc_gain(sca, H);
%
%
% % Now compute input impedance at Gate
% Zin_gate = sca.imp_looking_from_node(sca, 'G');





%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% 5T OTA
% sca = SCA_new();
%
%
%
% % Symbolic role parameters (one symbol per role)
% syms gm_in ro_in cgs_in cgd_in
% syms gm_ld ro_ld cgs_ld cgd_ld
% syms gm ro cgs cgd
% syms gm_tail ro_tail cgs_tail cgd_tail
% syms s Vin CL
%
%
% % Define differential input (this will create +Vin/2 and -Vin/2 vsrcs)
% % sca = sca.define_io(sca, {'in_p','in_n'}, {'out_p','out_n'}, 'Vin');
%
%
% sca = sca.define_io(sca, {'in_p','in_n'}, {'out_p'}, 'Vin');
%
%
%
% % % Define vdd input (PSR Analysis)
% % sca = sca.define_io(sca, {'vdd'}, {'out_p'}, 'Vin');
%
% mos = struct('gm',gm,'gmb',0,'ro',ro,'cgs',cgs,'cgd',0,'cgb',0,'cdb',0,'csb',0);
%
%
% % Canonical pins: {D, G, S, B}
% mos_in = struct('gm',gm_in,'gmb',0,'ro',ro,'cgs',cgs_in,'cgd',0,'cgb',0,'cdb',0,'csb',0);
%
% sca = sca.add_mos(sca,'M1', {'out_n' 'in_p','tail','0'}, mos_in);   % left input NMOS
% sca = sca.add_mos(sca,'M2', {'out_p','in_n','tail','0'}, mos_in);   % right input NMOS (matched)
%
%
% mos_ld = struct('gm',gm_ld,'gmb',0,'ro',ro,'cgs',cgs_ld,'cgd',0,'cgb',0,'cdb',0,'csb',0);
% sca = sca.add_mos(sca,'M3', {'out_n','out_n','vdd','vdd'}, mos_ld);
% sca = sca.add_mos(sca,'M4', {'out_p','out_n','vdd','vdd'}, mos_ld);
%
% % sca = sca.add_vsrc(sca,'Vbias_Load',{'mirror','vdd'}, 0);
% % sca = sca.add_vsrc(sca,'Vbias_inn',{'in_n','0'}, 0);
% % sca = sca.add_vsrc(sca,'Vbias_inp',{'in_p','0'}, 0);
% % sca = sca.add_vsrc(sca,'Vbias_tail',{'tail_bias','0'}, 0);
%
%
% % % TAIL DEVICE (M5)
% % mos_tail = struct('gm',gm_tail,'gmb',0,'ro',ro_tail,'cgs',cgs_tail,'cgd',0,'cgb',0,'cdb',0,'csb',0);
%
% % sca = sca.add_mos(sca,'M5', {'tail','tail_bias','0','0'}, mos_tail);
%
% % sca = sca.add_isrc(sca,'Ibias',{'tail','0'}, 0);
% sca = sca.add_res(sca,'RSS',{'tail','0'},sym('RSS'));
%
%
% sca = sca.define_node_role(sca,'vdd','supply');
% sca = sca.define_node_role(sca,'tail_bias','bias');
%
%
%
% % sca = sca.add_vsrc(sca,'VDD',{'vdd','0'}, 0);
%
%
% sca = sca.add_cap(sca,'C1',{'out_p','0'},CL);
% sca = sca.add_cap(sca,'C2',{'out_n','0'},CL);
%
%
%
% % Get transfer function
% H = sca.get_tf(sca);
%
% % Get DC gain
% Hdc = sca.get_dc_gain(sca, H);
%
% % % % Get common mode gain for differential input circuits
% H_cm = get_cm_gain(sca);
%
% % % % Get DC common mode gain
% Hdc = sca.get_dc_gain(sca, H_cm);
%
% % Calculate resistance looking into node
% Rin_out_p = sca.res_looking_from_node(sca, 'tail');
%
% % Calculate resistance looking into node
% Rin_out_n = sca.res_looking_from_node(sca, 'out_n');
%
% % % Simplify symbolic transfer functions using typical values
% % typical_values = struct('gm_ro', 50);
% % level = 10;
% % [H_simplified, info] = typicalValueSimplify(Hdc, typical_values, level);
%
% % % Calculate capacitance looking into node
% % Cap_1 = sca.cap_looking_from_node(sca, 'tail');
% %
% % % % Calculate capacitance looking into node
% % Cap_2 = sca.cap_looking_from_node(sca, 'out_n');
%
% show_netlist(sca);
%
% valid = sca.validate_netlist(sca);




%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Spectre Netlisting "" Two Stage Miller ""
% 
% sca = SCA_new();
% 
% % sca = sca.set_debug(sca, true);
% 
% % Call your parser
% sca = add_from_netlist_file(sca, 'Two_Stage_Miller_Netlist.txt');
% 
% sca = sca.update_component_param(sca, 'M6', 'ro', 'ro/2');
% 
% % sca = sca.update_component_param(sca, 'M4', 'ro', 'ro_4');
% % sca = sca.update_component_param(sca, 'M4', 'gm', 'gm_4');
% % sca = sca.update_component_param(sca, 'M2', 'gm', 'gm_2');
% % sca = sca.update_component_param(sca, 'M2', 'ro', 'ro_2');
% 
% sca = sca.define_io(sca, {'VINP', 'VINN'}, {'VOUT'}, 'Vin');
% 
% % res_short = sca.add_res(sca,'R_SC',{'VOUT','VINN'},0);
% 
% sca = sca.define_node_role(sca,'VDD','supply');
% sca = sca.define_node_role(sca,'GND','ground');
% 
% % Get transfer function
% [H, sca] = sca.get_tf(sca);
% 
% % Get DC gain
% Hdc = sca.get_dc_gain(sca, H);
% 
% % % Get common mode gain for differential input circuits
% % H_cm = get_cm_gain(sca);
% %
% % % Get DC common mode gain
% % Hdc = sca.get_dc_gain(sca, H_cm);
% 
% 
% % Calculate resistance looking into node
% Rin_Vout = sca.res_looking_from_node(sca, 'VOUT');
% 
% % Rin_out_FS = sca.res_looking_from_node(sca, 'net017');
% 
% Cin_Vout1 = sca.imp_looking_from_node(sca, 'net017', true);
% 
% Cin_Vout2 = sca.imp_looking_from_node(sca, 'VOUT', true);
% 
% % Get poles and zeros
% [poles, zeros] = sca.get_poles_zeros(sca, H);
% 
% 
% sca.show_netlist(sca);
% 
% valid = sca.validate_netlist(sca);



%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Spectre Netlisting "" Folded Cascode ""

% sca = SCA_new();
% 
% % sca = sca.set_debug(sca, true);
% 
% % Call your parser
% sca = add_from_netlist_file(sca, 'Folded_Cascode_Netlist.txt');
% 
% % Override MOSFET parameters
% sca = update_component_param(sca, 'M10', 'ro', 'ro/2');
% sca = update_component_param(sca, 'M12', 'ro', 'ro/2');
% 
% 
% sca = sca.define_io(sca, {'VINP', 'VINN'}, {'VOUTP', 'VOUTN'}, 'Vin');
% 
% 
% sca = sca.define_node_role(sca,'AVDD','supply');
% sca = sca.define_node_role(sca,'AGND','ground');
% sca = sca.define_node_role(sca,'VCASCP','bias');
% sca = sca.define_node_role(sca,'VCTRLP','bias');
% sca = sca.define_node_role(sca,'IB','bias');
% sca = sca.define_node_role(sca,'VBN','bias');
% sca = sca.define_node_role(sca,'VCASCN','bias');
% 
% 
% sca = sca.add_cap(sca,'C1',{'VOUTP','0'},sym('CL'));
% sca = sca.add_cap(sca,'C1',{'VOUTN','0'},sym('CL'));
% 
% % Get transfer function
% H = sca.get_tf(sca);
% 
% % Get DC gain
% Hdc = sca.get_dc_gain(sca, H);
% 
% 
% 
% % Get poles and zeros
% [poles, zeros] = sca.get_poles_zeros(sca, H);
% 
% 
% % Calculate resistance looking into node
% Rin_out_p = sca.res_looking_from_node(sca, 'VOUTP');
% 
% show_netlist(sca);
% 
% valid = sca.validate_netlist(sca);
% 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Cascode
% sca = SCA_new();
%
%
% % Symbolic role parameters
% syms gm_in ro_in cgs_in cgd_in
% syms gm_casc_p ro_casc_p cgs_casc_p cgd_casc_p
% syms gm_casc_n ro_casc_n cgs_casc_n cgd_casc_n
% syms gm_tail ro_tail cgs_tail cgd_tail
% syms s Vin CL RSS
%
%
% % Define input and output
% sca = sca.define_io(sca, {'in'}, {'out'}, 'Vin');
%
%
%
% % Define vdd input (PSR Analysis)
% % sca = sca.define_io(sca, {'vdd'}, {'out'}, 'Vin');
%
%
%
% mos = struct('gm',sym('gm'),'gmb',0,'ro',sym('ro'),'cgs',sym('cgs'),'cgd',0,'cgb',0,'cdb',0,'csb',0);
%
% mos_tail = struct('gm',gm_tail,'gmb',0,'ro',ro_tail,'cgs',cgs_tail,'cgd',cgd_tail,'cgb',0,'cdb',0,'csb',0);
%
% mos_casc_p = struct('gm',gm_casc_p,'gmb',0,'ro',ro_casc_p,'cgs',cgs_casc_p,'cgd',cgd_casc_p,'cgb',0,'cdb',0,'csb',0);
%
% mos_casc_n = struct('gm',gm_casc_n,'gmb',0,'ro',ro_casc_n,'cgs',cgs_casc_n,'cgd',cgd_casc_n,'cgb',0,'cdb',0,'csb',0);
%
% mos_in = struct('gm',gm_in,'gmb',sym('gmb_in'),'ro',ro_in,'cgs',cgs_in,'cgd',cgd_in,'cgb',0,'cdb',0,'csb',0);
%
%
% % Canonical pins: {D, G, S, B}
% sca = sca.add_mos(sca,'M1', {'node_p','0','vdd','vdd'}, mos);
%
% sca = sca.add_mos(sca,'M2', {'out','0','node_p','vdd'}, mos);
%
% sca = sca.add_mos(sca,'M3', {'out','Vbias_casc_n','node_n','0'}, mos);
%
% sca = sca.add_mos(sca,'M4', {'node_n','in','0','0'}, mos);
%
% % sca = sca.add_vsrc(sca,'Vbias_Load',{'mirror','vdd'}, 0);
% % sca = sca.add_vsrc(sca,'Vbias_tail',{'vdd', 'Vbias_tail'}, 0);
% % sca = sca.add_vsrc(sca,'Vbias_casc_p',{'vdd', 'Vbias_casc_p'}, 0);
% % sca = sca.add_vsrc(sca,'Vbias_casc_n',{'Vbias_casc_n','0'}, 0);
% % sca = sca.add_vsrc(sca,'VDD',{'vdd','0'}, 0);
%
% % Mark vdd as supply (AC-grounded unless used as IO)
% sca = sca.define_node_role(sca,'vdd','supply');
% sca = sca.define_node_role(sca,'Vbias_casc_n','bias');
%
%
% % sca = sca.add_isrc(sca,'Ibias',{'out','vdd'}, 0);
% % sca = sca.add_res(sca,'RSS',{'out','vdd'},RSS);
%
%
% sca = sca.add_cap(sca,'C1',{'out','0'},CL);
%
%
%
% % Get transfer function
% H = sca.get_tf(sca);
%
%
% % % Get poles and zeros
% % [poles, zeros] = sca.get_poles_zeros(sca, H);
%
% % Get DC gain
% Hdc = sca.get_dc_gain(sca, H);
%
% % % Calculate resistance looking into node
% % Rin_out_p = sca.res_looking_from_node(sca, 'node_n');
%
% % Calculate resistance looking into node
% % Rin_out_n = sca.res_looking_from_node(sca, 'out');
% Rin_out_n = sca.res_looking_from_node(sca, 'out');
% Rin_n = sca.res_looking_from_node(sca, 'node_n');
%
% % % Calculate capacitance looking into node
% % Cap_1 = sca.cap_looking_from_node(sca, 'out');
% %
% % % % Calculate capacitance looking into node
% % Cap_2 = sca.cap_looking_from_node(sca, 'node_n');
%
%
%
% % typical_vals = struct('gm_ro', 50);
% % [H_simp, info] = typicalValueSimplify(H, typical_vals, 10);
%
%
%
% show_netlist(sca);
%
% valid = sca.validate_netlist(sca);

%%%%%%%%%%%%%%%%%%%%%%%


% %% Spectre Netlisting "" BGR ""
%
% sca = SCA_new();
%
% sca = sca.set_debug(sca, true);
%
% % Call your parser
% sca = add_from_netlist_file(sca, 'BGR.txt');
%
% % sca = update_component_param(sca, 'M6', 'ro', 'ro/2');
%
% % sca = update_component_param(sca, 'M110', 'ro', 'ro/2');
% % sca = update_component_param(sca, 'M111', 'ro', 'ro/2');
% % sca = update_component_param(sca, 'M8', 'gm', 'gm_in');
% % sca = update_component_param(sca, 'M9', 'gm', 'gm_in');
%
% sca = sca.define_io(sca, {'vbe1', 'vbe2'}, {'vout_err'}, 'Vin');
%
% % res_short = sca.add_res(sca,'R_SC',{'VOUT','VINN'},0);
%
% sca = sca.define_node_role(sca,'SVDD','supply');
% sca = sca.define_node_role(sca,'vss','ground');
% sca = sca.define_node_role(sca,'vnmos','bias');
% sca = sca.define_node_role(sca,'net297','bias');
%
% % Get transfer function
% H = sca.get_tf(sca, true);
%
% % Get DC gain
% Hdc = sca.get_dc_gain(sca, H);
%
% % % % Get common mode gain for differential input circuits
% % % H_cm = get_cm_gain(sca);
% % %
% % % % Get DC common mode gain
% % % Hdc = sca.get_dc_gain(sca, H_cm);
% %
% %
% % % Calculate resistance looking into node
% % Rin_Vout = sca.res_looking_from_node(sca, 'VOUT');
% %
% % % Rin_out_FS = sca.res_looking_from_node(sca, 'net017');
% %
% % Cin_Vout1 = sca.cap_looking_from_node(sca, 'net017');
% %
% % Cin_Vout2 = sca.cap_looking_from_node(sca, 'VOUT');
% %
% % Get poles and zeros
% [poles, zeros] = sca.get_poles_zeros(sca, H);
%
%
% show_netlist(sca);
%
% valid = sca.validate_netlist(sca);













% Stop timer and record runtime
runtime = toc;

% Display runtime value
fprintf('Run Time = %.3f s\n\n',runtime);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  *********************************                      %
%                  *   __  __   ______   _____     *                      %
%                  *  |  \/  | |  ____| |  __ \    *                      %
%                  *  | \  / | | |____  | |__| |   *                      %
%                  *  | |\/| | |  ____| |  __  |   *                      %
%                  *  | |  | | | |      | |  | |   *                      %
%                  *  |_|  |_| |_|      |_|  |_|   *                      %
%                  *                               *                      %
%                  *********************************                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Symbolic Circuit Analyzer (SCA)
% - Builds netlist (struct)
% - Stamps passive and MOS small-signal pi model into symbolic G and C
% - Supports independent voltage and current sources.
% - Extracts symbolic transfer H(s) = Vout/Vin
% - Handles supply and bias nodes (auto-merged to ground unless defined as IO)
% - Added error checking, validation, and utility functions
%


function sca = SCA_new()
% Initialize an empty SCA struct
sca.nodes = {'0'}; % ground node is '0'
sca.node_map = containers.Map({'0'},{1}); % name->index (1-based)
sca.components = {}; % list of component structs
sca.vsrc_count = 0;
sca.Nnodes = 1; % count including ground
sca.s = sym('s');
sca.io_nodes = {}; % track IO nodes
sca.roles = containers.Map(); % node roles (supply, etc.)
sca.debug_mode = false; % for debugging output
sca = classify(sca); % attach methods

% analysis matrices (initialize empty)
sca.matrices = {};

end


function sca = classify(sca)
% Attach functions
sca.add_res = @add_res;
sca.add_cap = @add_cap;
sca.add_vsrc = @add_vsrc;
sca.add_vccs = @add_vccs;
sca.add_isrc = @add_isrc;
sca.add_mos = @add_mos;
sca.define_io = @define_io;
sca.define_node_role = @define_node_role;
sca.merge_nodes = @merge_nodes;
sca.nodemap = @nodemap;
sca.get_tf = @get_tf;
sca.get_dc_gain = @get_dc_gain;
sca.get_poles_zeros = @get_poles_zeros;
sca.get_cm_gain = @get_cm_gain;
sca.imp_looking_from_node = @imp_looking_from_node;
sca.res_looking_from_node = @res_looking_from_node;
sca.cap_looking_from_node = @cap_looking_from_node;
sca.add_from_netlist = @add_from_netlist;
sca.add_from_netlist_file = @add_from_netlist_file;
sca.typicalValueSimplify = @typicalValueSimplify;
sca.update_component_param = @update_component_param;
sca.show_netlist = @show_netlist;
sca.validate_netlist = @validate_netlist;
sca.set_debug = @set_debug;
end


function sca = add_node_if_missing(sca, name)
if isempty(name) || strcmp(name,'')
    error('SCA:EmptyNodeName', 'Empty node name not allowed');
end
if ~isKey(sca.node_map, name)
    sca.Nnodes = sca.Nnodes + 1;
    sca.nodes{end+1} = name;
    sca.node_map(name) = sca.Nnodes;
    if sca.debug_mode
        fprintf('Added node: %s (index %d)\n', name, sca.Nnodes);
    end
end
end


function sca = add_res(sca, name, pins, Rval)
if numel(pins) ~= 2
    error('SCA:InvalidPins', 'Resistor requires exactly 2 pins');
end

% Ensure nodes exist in the circuit
sca = add_node_if_missing(sca, pins{1});
sca = add_node_if_missing(sca, pins{2});

% Define component structure
comp.type = 'R';
comp.name = name;
comp.pins = pins;
comp.params.R = sym(Rval);
sca.components{end+1} = comp;
if sca.debug_mode
    fprintf('Added resistor %s: %s-%s, R=%s\n', name, pins{1}, pins{2}, char(comp.params.R));
end
end


function sca = add_cap(sca, name, pins, Cval)
if numel(pins) ~= 2
    error('SCA:InvalidPins', 'Capacitor requires exactly 2 pins');
end

% Ensure nodes exist in the circuit
sca = add_node_if_missing(sca, pins{1});
sca = add_node_if_missing(sca, pins{2});

% Define component structure
comp.type = 'C';
comp.name = name;
comp.pins = pins;
comp.params.C = sym(Cval);
sca.components{end+1} = comp;
if sca.debug_mode
    fprintf('Added capacitor %s: %s-%s, C=%s\n', name, pins{1}, pins{2}, char(comp.params.C));
end
end

function sca = add_ind(sca, name, pins, Lval)
if numel(pins) ~= 2
    error('SCA:InvalidPins', 'Inductor requires exactly 2 pins');
end

% Ensure nodes exist in the circuit
sca = add_node_if_missing(sca, pins{1});
sca = add_node_if_missing(sca, pins{2});

% Define component structure
comp.type = 'L';
comp.name = name;
comp.pins = pins;
comp.params.L = sym(Lval);

% Add to circuit
sca.components{end+1} = comp;

% Debug print
if sca.debug_mode
    fprintf('Added inductor %s: %s-%s, L=%s\n', name, pins{1}, pins{2}, char(comp.params.L));
end
end


function sca = add_vsrc(sca, name, pins, amp_sym)
if numel(pins) ~= 2
    error('SCA:InvalidPins', 'Voltage source requires exactly 2 pins');
end
sca = add_node_if_missing(sca, pins{1});
sca = add_node_if_missing(sca, pins{2});
sca.vsrc_count = sca.vsrc_count + 1;
comp.type = 'V';
comp.name = name;
comp.pins = pins;
comp.params.amp = sym(amp_sym);
comp.params.index = sca.vsrc_count;
sca.components{end+1} = comp;
if sca.debug_mode
    fprintf('Added voltage source %s: %s-%s, V=%s\n', name, pins{1}, pins{2}, char(comp.params.amp));
end
end


function sca = add_isrc(sca, name, pins, amp_sym)
% Ideal current source definition
% pins = {positive_node, negative_node}
comp.type = 'I';
comp.name = name;
comp.pins = pins;
comp.params.amp = amp_sym;  % symbolic or numeric amplitude
sca.components{end+1} = comp;
if sca.debug_mode
    fprintf('Added current source %s: %s-%s, I=%s\n', name, pins{1}, pins{2}, char(comp.params.amp));
end
end

function sca = add_vccs(sca, name, out_nodes, ctrl_nodes, value)
% Add a voltage-controlled current source
% out_nodes = {nplus, nminus} current flows from nplus → nminus
% ctrl_nodes = {cplus, cminus} controlling voltage V(cplus,cminus)

nplus = nodemap(sca, out_nodes{1});
nminus = nodemap(sca, out_nodes{2});
cplus = nodemap(sca, ctrl_nodes{1});
cminus = nodemap(sca, ctrl_nodes{2});

g = value;

% Stamping into G matrix:
% I(out_nodes) = g * (V(cplus) - V(cminus))

if nplus ~= 0
    if cplus ~= 0, sca.G(nplus,cplus) = sca.G(nplus,cplus) + g; end
    if cminus ~= 0, sca.G(nplus,cminus) = sca.G(nplus,cminus) - g; end
end
if nminus ~= 0
    if cplus ~= 0, sca.G(nminus,cplus) = sca.G(nminus,cplus) - g; end
    if cminus ~= 0, sca.G(nminus,cminus) = sca.G(nminus,cminus) + g; end
end
end


function sca = add_mos(sca, name, pins, params)
% pins = {D,G,S,B} or {D,G,S} (B defaults to '0')
if numel(pins) < 3
    error('SCA:InvalidPins', 'MOSFET requires at least D,G,S pins');
end

% Default bulk to ground if not specified
if numel(pins) == 3
    pins{4} = '0';
end

for k=1:numel(pins)
    sca = add_node_if_missing(sca, pins{k});
end

% Default MOSFET parameters
default_params = struct('gm',0, 'gmb',0, 'ro',0, 'cgs',0, 'cgd',0, 'cgb',0, 'cdb',0, 'csb',0);
param_names = fieldnames(default_params);

for i=1:length(param_names)
    pname = param_names{i};
    if ~isfield(params, pname)
        params.(pname) = sym(0);
    else
        params.(pname) = sym(params.(pname));
    end
end

comp.type = 'M';
comp.name = name;
comp.pins = pins;
comp.params = params;
sca.components{end+1} = comp;

if sca.debug_mode
    fprintf('Added MOSFET %s: D=%s, G=%s, S=%s, B=%s\n', name, pins{1}, pins{2}, pins{3}, pins{4});
end
end


function sca = define_io(sca, in_nodes, out_nodes, in_sym)
% define_io - Define input and output nodes for the circuit
%
% Automatically determines if inputs/outputs are single-ended or differential
% based on the number of nodes provided:
% - 1 node  = single-ended
% - 2 nodes = differential
%
% Inputs:
%   sca       - SCA circuit object
%   in_nodes  - Cell array of input node names (1 for SE, 2 for diff)
%   out_nodes - Cell array of output node names (1 for SE, 2 for diff)
%   in_sym    - Symbolic variable for input voltage
%
% Output:
%   sca       - Updated SCA circuit object

% Validate inputs
if ~iscell(in_nodes)
    if ischar(in_nodes) || isstring(in_nodes)
        in_nodes = {char(in_nodes)};
    else
        error('SCA:InvalidInput', 'in_nodes must be a cell array of node names');
    end
end

if ~iscell(out_nodes)
    if ischar(out_nodes) || isstring(out_nodes)
        out_nodes = {char(out_nodes)};
    else
        error('SCA:InvalidInput', 'out_nodes must be a cell array of node names');
    end
end

% Determine input type based on number of nodes
num_in_nodes = numel(in_nodes);
if num_in_nodes == 1
    in_type = 'single-ended';
elseif num_in_nodes == 2
    in_type = 'differential';
else
    error('SCA:InvalidInputNodes', ...
        'Input must have 1 node (single-ended) or 2 nodes (differential), rot %d nodes', ...
        num_in_nodes);
end

% Determine output type based on number of nodes
num_out_nodes = numel(out_nodes);
if num_out_nodes == 1
    out_type = 'single-ended';
elseif num_out_nodes == 2
    out_type = 'differential';
else
    error('SCA:InvalidOutputNodes', ...
        'Output must have 1 node (single-ended) or 2 nodes (differential), rot %d nodes', ...
        num_out_nodes);
end

% Validate that all nodes exist or can be created
all_nodes = [in_nodes, out_nodes];
for i = 1:length(all_nodes)
    if ~any(strcmp(sca.nodes, all_nodes{i}))
        if sca.debug_mode
            fprintf('Adding missing I/O node: %s\n', all_nodes{i});
        end
        sca = add_node_if_missing(sca, all_nodes{i});
    end
end

% Store I/O configuration
sca.in_nodes = in_nodes;
sca.out_nodes = out_nodes;
sca.in_type = in_type;
sca.out_type = out_type;
sca.in_sym = sym(in_sym);

% Update I/O nodes list (avoid duplicates)
new_io_nodes = [in_nodes, out_nodes];
for i = 1:length(new_io_nodes)
    if ~any(strcmp(sca.io_nodes, new_io_nodes{i}))
        sca.io_nodes{end+1} = new_io_nodes{i};
    end
end

% Add input voltage sources based on determined type
if strcmp(in_type, 'differential')
    % Differential input: Vin_p = +Vin/2, Vin_n = -Vin/2
    sca = add_vsrc(sca, 'Vin_p', {in_nodes{1}, '0'}, sca.in_sym/2);
    sca = add_vsrc(sca, 'Vin_n', {in_nodes{2}, '0'}, -sca.in_sym/2);
    
    if sca.debug_mode
        fprintf('Added differential input sources:\n');
        fprintf('  Vin_p: %s to ground = +%s/2\n', in_nodes{1}, char(sca.in_sym));
        fprintf('  Vin_n: %s to ground = -%s/2\n', in_nodes{2}, char(sca.in_sym));
    end
else
    % Single-ended input
    sca = add_vsrc(sca, 'Vin', {in_nodes{1}, '0'}, sca.in_sym);
    
    if sca.debug_mode
        fprintf('Added single-ended input source:\n');
        fprintf('  Vin: %s to ground = %s\n', in_nodes{1}, char(sca.in_sym));
    end
end

% Set node roles for better circuit understanding
for i = 1:length(in_nodes)
    sca = define_node_role(sca, in_nodes{i}, 'input');
end
for i = 1:length(out_nodes)
    sca = define_node_role(sca, out_nodes{i}, 'output');
end

% Debug output
if sca.debug_mode
    fprintf('Defined I/O configuration:\n');
    fprintf('  Input:  %s (%s)\n', in_type, strjoin(in_nodes, ', '));
    fprintf('  Output: %s (%s)\n', out_type, strjoin(out_nodes, ', '));
    fprintf('  Input symbol: %s\n', char(sca.in_sym));
end
end


function sca = define_node_role(sca, node, role)
% define_node_role - Assign a role to a node
%
% Inputs:
%   sca  - SCA circuit object
%   node - Node name (string/char)
%   role - Node role (string/char)
%
% Output:
%   sca  - Updated SCA circuit object

valid_roles = {'supply', 'ground', 'input', 'output', 'internal', 'bias'};

if ~ismember(role, valid_roles)
    warning('SCA:UnknownRole', 'Unknown node role: %s', role);
end

% Ensure node exists
sca = add_node_if_missing(sca, node);

% Initialize roles container if it doesn't exist
if ~isfield(sca, 'roles')
    sca.roles = containers.Map();
end

% Set the role
sca.roles(node) = role;

% Handle special roles
switch role
    case 'ground'
        % Force merge to global ground
        sca = sca.merge_nodes(sca, node, '0');
        
    case {'supply', 'bias'}
        % By default these should be AC grounded unless explicitly I/O
        if ~ismember(node, sca.io_nodes)
            if sca.debug_mode
                fprintf('Merging %s node %s to ground (AC reference)\n', role, node);
            end
            sca = sca.merge_nodes(sca, node, '0');
        end
        
    case 'input'
        % Inputs are left floating unless explicitly tied elsewhere
        if sca.debug_mode
            fprintf('Node %s set as INPUT\n', node);
        end
        
    case 'output'
        % Outputs are observable nodes
        if sca.debug_mode
            fprintf('Node %s set as OUTPUT\n', node);
        end
        
    case 'internal'
        % Internal nodes are untouched
        if sca.debug_mode
            fprintf('Node %s set as INTERNAL\n', node);
        end
end

if sca.debug_mode
    fprintf('Assigned role "%s" to node "%s"\n', role, node);
end
end


function sca = merge_nodes(sca, node1, node2)
% Rewire all references of node1 to node2
if strcmp(node1, node2)
    return; % nothing to do
end

if sca.debug_mode
    fprintf('Merging node %s into %s\n', node1, node2);
end

for k=1:numel(sca.components)
    pins = sca.components{k}.pins;
    for p=1:numel(pins)
        if strcmp(pins{p}, node1)
            pins{p} = node2;
        end
    end
    sca.components{k}.pins = pins;
end

% Remove node1 from map and list
if isKey(sca.node_map, node1)
    remove(sca.node_map, node1);
end
sca.nodes(strcmp(sca.nodes, node1)) = [];

% Update node indices
sca.node_map = containers.Map();
for i=1:length(sca.nodes)
    sca.node_map(sca.nodes{i}) = i;
end
sca.Nnodes = length(sca.nodes);
end


function show_netlist(sca)
fprintf('\n=== SCA Netlist ===\n');
fprintf('Nodes (%d): ', sca.Nnodes);
for i=1:length(sca.nodes)
    fprintf('%s(%d) ', sca.nodes{i}, i);
end
fprintf('\n');

if ~isempty(keys(sca.roles))
    fprintf('Node roles:\n');
    role_keys = keys(sca.roles);
    for i=1:length(role_keys)
        fprintf('  %s: %s\n', role_keys{i}, sca.roles(role_keys{i}));
    end
end

fprintf('Components (%d):\n', length(sca.components));
for k=1:numel(sca.components)
    c = sca.components{k};
    fprintf('  %s %s: ', c.type, c.name);
    fprintf('pins=[%s] ', strjoin(c.pins, ','));
    
    % Show key parameters
    switch c.type
        case 'R'
            fprintf('R=%s', char(c.params.R));
        case 'C'
            fprintf('C=%s', char(c.params.C));
        case 'V'
            fprintf('V=%s', char(c.params.amp));
        case 'M'
            fprintf('gm=%s ro=%s', char(c.params.gm), char(c.params.ro));
        case 'isrc'
            pos = nidx(comp.pins{1});
            neg = nidx(comp.pins{2});
            I = comp.params.amp;
            if pos ~= 0, B(pos) = B(pos) - I; end   % current leaving pos node
            if neg ~= 0, B(neg) = B(neg) + I; end   % current entering neg node
    end
    fprintf('\n');
end

if isfield(sca, 'in_nodes')
    fprintf('Input: %s (%s)\n', strjoin(sca.in_nodes, ','), sca.in_type);
    fprintf('Output: %s (%s)\n', strjoin(sca.out_nodes, ','), sca.out_type);
end
fprintf('==================\n\n');
end


function valid = validate_netlist(sca)
% Basic netlist validation
valid = true;

% Check if I/O is defined
if ~isfield(sca, 'in_nodes') || ~isfield(sca, 'out_nodes')
    fprintf('Warning: I/O not defined\n');
    valid = false;
end

% Check for floating nodes (nodes with only one connection)
node_connections = containers.Map();
for node = sca.nodes
    node_connections(node{1}) = 0;
end

for k=1:numel(sca.components)
    pins = sca.components{k}.pins;
    for pin = pins
        if isKey(node_connections, pin{1})
            node_connections(pin{1}) = node_connections(pin{1}) + 1;
        end
    end
end

floating_nodes = {};
node_keys = keys(node_connections);
for i=1:length(node_keys)
    node = node_keys{i};
    if node_connections(node) <= 1 && ~strcmp(node, '0')
        floating_nodes{end+1} = node;
    end
end

if ~isempty(floating_nodes)
    fprintf('Warning: Floating nodes detected: %s\n', strjoin(floating_nodes, ', '));
end

fprintf('Netlist validation %s\n', ternary(valid, 'passed', 'failed'));
end


function sca = set_debug(sca, debug_flag)
sca.debug_mode = debug_flag;
fprintf('Debug mode: %s\n', ternary(debug_flag, 'ON', 'OFF'));
end


function idx = nodemap(sca ,name)
% Node mapping function with AC-ground handling for special roles

% Treat '0' as ground
if strcmp(name,'0')
    idx = 0;
    return;
end

% If roles container exists and this node is marked as supply/bias
% and the node is NOT an I/O node, treat it as AC ground (index 0).
try
    if isfield(sca,'roles') && isKey(sca.roles, name)
        role = sca.roles(name);
        if (strcmp(role,'supply') || strcmp(role,'bias'))
            % If io_nodes exists and node is part of I/O, do not ground it
            if isfield(sca,'io_nodes') && ~isempty(sca.io_nodes) && any(strcmp(sca.io_nodes, name))
                % leave as node (do nothing here)
            else
                % AC-ground this node
                if sca.debug_mode
                    fprintf('AC-grounding node "%s" (role=%s) for small-signal analysis\n', name, char(role));
                end
                idx = 0;
                return;
            end
        end
    end
catch
    % If roles map access fails for any reason, fall back to normal behavior
end

% Normal mapping via node_map
if ~isKey(sca.node_map, name)
    error('Node %s not found in node_map', name);
end
n = sca.node_map(name);
if n == 1
    idx = 0;
else
    idx = n - 1;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Stamping functions
function G = stampR(n1, n2, R, G)
if R == 0 || R == sym(0)
    return; % Skip zero resistance (short circuit)
end
g = 1/R;
if ~isequal(n1, 0)
    G(n1,n1) = G(n1,n1) + g;
end
if ~isequal(n2, 0)
    G(n2,n2) = G(n2,n2) + g;
end
if ~isequal(n1, 0) && ~isequal(n2, 0)
    G(n1,n2) = G(n1,n2) - g;
    G(n2,n1) = G(n2,n1) - g;
end
end

function C = stampC(n1, n2, Cval, C)
if isequal(Cval, 0) || isequal(Cval, sym(0))
    return; % Skip zero capacitance
end
if ~isequal(n1, 0)
    C(n1,n1) = C(n1,n1) + Cval;
end
if ~isequal(n2, 0)
    C(n2,n2) = C(n2,n2) + Cval;
end
if ~isequal(n1, 0) && ~isequal(n2, 0)
    C(n1,n2) = C(n1,n2) - Cval;
    C(n2,n1) = C(n2,n1) - Cval;
end
end


function [G, B] = stampV(nplus, nminus, amp, N, vs_idx, G, B)
row = N + vs_idx;
if nplus
    G(nplus,row) = G(nplus,row) + 1;
    G(row,nplus) = G(row,nplus) + 1;
end
if nminus
    G(nminus,row) = G(nminus,row) - 1;
    G(row,nminus) = G(row,nminus) - 1;
end
B(row) = amp;

end


function G = stampVCCS(p, n, pc, nc, g, G)

% out_nodes = {p, n} current flows from p → n
% ctrl_nodes = {pc, nc} controlling voltage V(pc,nc)

% VCCS stamp:
% i = g * (v_pc - v_nc)
if p && pc
    G(p,pc) = G(p,pc) + g;
end
if p && nc
    G(p,nc) = G(p,nc) - g;
end
if n && pc
    G(n,pc) = G(n,pc) - g;
end
if n && nc
    G(n,nc) = G(n,nc) + g;
end

end


function [G, C] = stampMOS(Drain, Gate, Source, Bulk, par, G, C)

% Output resistance (ro between Drain and Source)
G = stampR(Drain,Source,par.ro,G);

% gm current: gm·vgs current from Drain to Source (controlled by Gate-Source)
G = stampVCCS(Drain,Source,Gate,Source,par.gm,G);

% gmb current: gmb·vbs current from Drain to Source (controlled by Bulk-Source)
G = stampVCCS(Drain,Source,Bulk,Source,par.gmb,G);

% Capacitors
C = stampC(Gate,Source,par.cgs,C);
C = stampC(Gate,Drain,par.cgd,C);
C = stampC(Drain,Bulk,par.cdb,C);
C = stampC(Gate,Bulk,par.cgb,C);
C = stampC(Source,Bulk,par.csb,C);

end


function B = stampI(nplus, nminus, amp, B)

try
    if nplus ~= 0, B(nplus) = B(nplus) - amp; end   % current leaving nplus node
    if nminus ~= 0, B(nminus) = B(nminus) + amp; end   % current entering nminus node
catch
    % ignore malformed isrc
end
end


function [G, C, B] = stamp_all_components(sca, N, G, C, B)

vs_idx = 0;
for k=1:numel(sca.components)
    c = sca.components{k};
    switch c.type
        case 'R'
            G = stampR(nodemap(sca, c.pins{1}), nodemap(sca, c.pins{2}), c.params.R, G);
        case 'C'
            C = stampC(nodemap(sca, c.pins{1}), nodemap(sca, c.pins{2}), c.params.C, C);
        case 'isrc'
            B = stampI(nodemap(sca, c.pins{1}), nodemap(sca, c.pins{2}), c.params.amp, B);
        case 'V'
            vs_idx = vs_idx + 1;
            [G, B] = stampV(nodemap(sca, c.pins{1}), nodemap(sca, c.pins{2}), c.params.amp, N, vs_idx, G, B);
        case 'VCCS'
            vs_idx = vs_idx + 1;
            [G] = stampVCCS(nodemap(sca, c.pins{1}), nodemap(sca, c.pins{2}), nodemap(sca, c.pins{3}), nodemap(sca, c.pins{4}), c.params.g, G);
        case 'M'
            pins = c.pins;
            while numel(pins) < 4
                pins{end+1} = '0';
            end
            [G, C] = stampMOS(nodemap(sca, pins{1}), nodemap(sca, pins{2}), nodemap(sca, pins{3}),...
                nodemap(sca, pins{4}), c.params, G, C);
            
    end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'


function [H, sca] = get_tf(sca , print_flag)

if nargin < 2
    print_flag = false;
end

if sca.debug_mode
    fprintf('Computing transfer function...\n');
    tic;
end

s = sca.s;

% Count voltage sources
Nv = 0;
for k=1:numel(sca.components)
    if strcmp(sca.components{k}.type,'V')
        Nv = Nv + 1;
    end
end

N = sca.Nnodes - 1; % excluding ground
M = N + Nv;
G = sym(zeros(M,M));
C = sym(zeros(M,M));
B = sym(zeros(M,1));



% Stamp all components
[G, C, B] = stamp_all_components(sca, N, G, C, B);


% Saving MNA matrices to be used if needed
MNA_matrices.G = G;
MNA_matrices.C = C;
MNA_matrices.B = B;
sca.matrices = MNA_matrices;


if sca.debug_mode
    fprintf('Matrix size: %dx%d\n', size(G,1), size(G,2));
end

% Solve for transfer function
Y = G + s*C;
X = Y \ B;


if sca.debug_mode
    fprintf('Rank(G+sC): %d (of %d)\n', rank(Y), size(Y,1));
end



% Extract output voltage
if strcmp(sca.out_type,'single-ended')
    idx = nodemap(sca, sca.out_nodes{1});
    Vout = (idx==0)*0 + (idx~=0)*X(idx);
else
    idxp = nodemap(sca, sca.out_nodes{1});
    idxn = nodemap(sca, sca.out_nodes{2});
    Vout = (idxp~=0)*X(idxp) - (idxn~=0)*X(idxn);
end

try
    H = simplify(Vout / sca.in_sym);
    
catch
    H = Vout / sca.in_sym;
end

if sca.debug_mode
    elapsed = toc;
    fprintf('Transfer function computed in %.3f seconds\n', elapsed);
end

fprintf('Transfer function H(s) computed successfully.\n');

if print_flag
    fprintf('H(s) =\n');
    pretty(H);
end

end


function Hdc = get_dc_gain(sca, Hs)
% Get DC gain of transfer function by setting s=0
if nargin < 2
    tf_temp = sca.get_tf(sca);
    Hs = tf_temp;
end
Hdc = simplify(subs(Hs, sym('s'), 0));
Hdc = simplifyFraction(Hdc);
fprintf('DC gain:\n');
pretty(Hdc);
end


function [poles, zeros] = get_poles_zeros(sca, Hs)
% Extract poles and zeros from transfer function
if nargin < 2
    Hs = sca.get_tf(sca);
end

try
    [num, den] = numden(Hs);
    
    % Find zeros (roots of numerator)
    zeros_coeffs = coeffs(num, sca.s);
    if length(zeros_coeffs) > 1
        zeros = solve(num == 0, sca.s);
    else
        zeros = [];
    end
    
    % Find poles (roots of denominator)
    poles_coeffs = coeffs(den, sca.s);
    if length(poles_coeffs) > 1
        poles = solve(den == 0, sca.s);
    else
        poles = [];
    end
    
    fprintf('Poles:\n');
    if ~isempty(poles)
        pretty(poles);
    else
        fprintf('No poles found (constant transfer function)\n\n');
    end
    
    fprintf('Zeros:\n');
    if ~isempty(zeros)
        pretty(zeros);
    else
        fprintf('No zeros found\n\n');
    end
    
catch ME
    fprintf('Could not extract poles/zeros: %s\n\n', ME.message);
    poles = [];
    zeros = [];
end
end


function Hcm = get_cm_gain(sca, print_flag)
% Get common mode gain for differential input circuits
% Both inputs are set to the same value (Vcm)

if nargin < 2
    print_flag = false;
end

if ~isfield(sca, 'in_type') || ~strcmp(sca.in_type, 'differential')
    error('SCA:NotDifferential', 'Common mode gain only applies to differential inputs');
end

if sca.debug_mode
    fprintf('Computing common mode gain...\n');
    tic;
end

% Create a temporary SCA with common mode sources
temp_sca = sca;
Vcm = sym('Vcm'); % Common mode input voltage

% Remove existing differential voltage sources and add common mode ones
temp_components = {};
for k = 1:numel(sca.components)
    comp = sca.components{k};
    % Skip the original differential input sources
    if strcmp(comp.type, 'V') && (strcmp(comp.name, 'Vin_p') || strcmp(comp.name, 'Vin_n'))
        continue;
    end
    temp_components{end+1} = comp;
end

% Add new common mode voltage sources (both inputs at same potential)
temp_sca.components = temp_components;
temp_sca = add_vsrc(temp_sca, 'Vcm_p', {sca.in_nodes{1}, '0'}, Vcm);
temp_sca = add_vsrc(temp_sca, 'Vcm_n', {sca.in_nodes{2}, '0'}, Vcm);

% Calculate transfer function with common mode input
temp_sca.in_sym = Vcm; % Update input symbol for transfer function calculation
Hcm = temp_sca.get_tf(temp_sca, print_flag);

if sca.debug_mode
    elapsed = toc;
    fprintf('Common mode gain computed in %.3f seconds\n', elapsed);
end

if print_flag
    fprintf('Common mode gain H_CM(s):\n');
    pretty(Hcm);
end

end


% Utility function
function result = ternary(condition, true_val, false_val)
if condition
    result = true_val;
else
    result = false_val;
end
end


%%%%%%%%%%%%%%%%%%%%%%%


function Zin = imp_looking_from_node(sca, node_name, print_flag)
% imp_looking_from_node - Symbolic input impedance at a specified node
%
% Zin(s) = Vtest / Itest  with 1A current source injected at node_name


if nargin < 3
    print_flag = false;
end

% --- Validate node ---
if ~ischar(node_name) && ~isstring(node_name)
    error('node_name must be a string or char array');
end
node_name = char(node_name);

if ~any(strcmp(sca.nodes, node_name))
    error('Node "%s" not found in circuit. Available nodes: %s', ...
        node_name, strjoin(sca.nodes, ', '));
end

node_idx = nodemap(sca, node_name);
if node_idx == 0
    error('Cannot calculate capacitance looking into ground node');
end

% --- Backup ---
original_components = sca.components;
s = sca.s;

% --- Zero out independent DC sources ---
for k = 1:length(sca.components)
    comp = sca.components{k};
    if isfield(comp,'type')
        t = upper(comp.type);
        if startsWith(t,'V')
            if isfield(comp,'params') && isfield(comp.params,'amp')
                sca.components{k}.params.amp = sym(0);
            end
        elseif startsWith(t,'I')
            if isfield(comp,'params') && isfield(comp.params,'amp')
                sca.components{k}.params.amp = sym(0);
            end
        end
    end
end

% --- Count voltage sources for MNA size ---
N = sca.Nnodes - 1;
Nv = 0;
for k = 1:numel(sca.components)
    c = sca.components{k};
    if isfield(c,'type') && (strcmpi(c.type,'V') || startsWith(upper(c.type),'V'))
        Nv = Nv + 1;
    end
end
M = N + Nv;

G = sym(zeros(M,M));
C = sym(zeros(M,M));
B = sym(zeros(M,1));

% Stamp all components
[G, C, B] = stamp_all_components(sca, N, G, C, B);

% --- Build AC MNA matrix ---
Y = G + s*C;

% --- Inject 1A test current ---
Itest = sym(zeros(size(Y,1),1));
Itest(node_idx) = 1;


% --- Solve voltages ---
V = Y \ Itest;

% --- Extract driving-point functions ---
Vin = V(node_idx);
Zin = simplify(Vin);   % master impedance

% Cancel common factors in a rational expression
[N,D] = numden(Zin);
Zin = simplifyFraction(N/D);  % cancel gcd between numerator & denominator

if sca.debug_mode
    fprintf('Impedance looking into "%s":\n', node_name);
    pretty(Zin);
end

if print_flag
    fprintf('Impedance looking into "%s":\n', node_name);
    pretty(Zin);
end
end

function Rin = res_looking_from_node(sca, node_name)
% res_looking_from_node - Equivalent resistance at node (Re{Zin(s)} @ s=0)

syms s
Zin = sca.imp_looking_from_node(sca, node_name);

% Extract DC resistance
% Rin = simplify(limit(real(Zin), s, 0));
Rin = simplify(limit(Zin, s, 0));

fprintf('Resistance looking into "%s":\n', node_name);
pretty(Rin);

end

function Cin = cap_looking_from_node(sca, node_name)
% cap_looking_from_node - Equivalent capacitance at node (dYin/ds @ s=0)

syms s
Zin = sca.imp_looking_from_node(sca, node_name);

% Admittance
Yin = simplify(1/Zin);

% Small-s expansion: Y(s) ≈ G0 + s*Ceq + ...
dYds = diff(Yin, s);
Cin = simplify(subs(dYds, s, 0));

%%%%%%%%%%%%%%%%%%
% Cin = partfrac(Yin, s);
% Cin = simplify(partfrac(Yin, s));


%%%%%%%%%%%%%%%%%%
% % Yin is symbolic admittance (function of s)
% w = sym('w','real');
% Yjw = simplify(subs(Yin, s, 1i*w));      % substitute s = j*w
% Cin = simplify(imag(Yjw)/w);% imag(Y)/w -> C at low freq
%
% % tidy the result: cancel common factors in a rational expression
% [N,D] = numden(Cin);
% Cin = simplifyFraction(N/D);  % cancel gcd between numerator & denominator


fprintf('Capacitance looking into "%s":\n', node_name);
pretty(Cin);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sca = add_from_netlist_file(sca, filename, ignore_sources)
% Read Cadence Spectre netlist from .txt file and build SCA circuit
%
% Usage: sca = add_from_netlist_file(sca, 'netlist.txt')
%
% Features:
% - Ignores MOSFET models and device sizing
% - Ignores voltage and current sources "optinal"
% - Uses 0V default for voltage sources (ignores dc values)
% - Extracts node connectivity
% - Assigns default small-signal parameters for MOSFETs

% Add optional parameter with default value "true"
if nargin < 3
    ignore_sources = true;
end

try
    % Read file content
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    % Read all lines
    netlistLines = {};
    lineCount = 0;
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)  % fgetl returns -1 at EOF
            lineCount = lineCount + 1;
            netlistLines{lineCount} = line;
        end
    end
    fclose(fid);
    
    fprintf('Successfully read %d lines from %s\n', lineCount, filename);
    
    % Process the netlist
    sca = add_from_netlist(sca, netlistLines, ignore_sources);
    
catch ME
    if exist('fid', 'var') && fid ~= -1
        fclose(fid);
    end
    rethrow(ME);
end


end

function sca = add_from_netlist(sca, netlistText, ignore_sources)
% Process netlist text (cell array of strings or single string)

% Add optional parameter with default value "true"
if nargin < 3
    ignore_sources = true;
end

if ischar(netlistText) || isstring(netlistText)
    % Single string, split by newlines
    if ischar(netlistText)
        netlistText = string(netlistText);
    end
    lines = splitlines(netlistText);
elseif iscell(netlistText)
    % Cell array of strings
    lines = string(netlistText);
else
    error('Unsupported netlist format');
end

% Remove empty entries
lines = lines(strlength(lines) > 0);

% Process each line
processedCount = 0;
for k = 1:numel(lines)
    line = strtrim(lines(k));
    
    % Skip empty lines and comments
    if strlength(line) == 0 || startsWith(line, '*') || startsWith(line, '//')
        continue;
    end
    
    % Skip simulator commands and subckt definitions
    if startsWith(line, '.') || startsWith(line, 'simulator') || ...
            contains(line, 'subckt') || contains(line, 'ends')
        continue;
    end
    
    try
        sca = add_from_netlist_line(sca, line, ignore_sources);
        processedCount = processedCount + 1;
    catch ME
        warning('Failed to process line %d: %s\nError: %s', k, char(line), ME.message);
    end
end

fprintf('Successfully processed %d circuit elements\n', processedCount);
end

function sca = add_from_netlist_line(sca, line, ignore_sources)
% Parse a single line of Spectre netlist and add to SCA

% Add optional parameter with default value "true"
if nargin < 3
    ignore_sources = true;
end


% Convert to char for easier string manipulation
line = char(line);

% Split line into tokens
tokens = strsplit(strtrim(line));
if length(tokens) < 2
    return;  % Skip invalid lines
end

name = tokens{1};           % Instance name
compType = upper(name(1));  % Component type from first character

% Extract nodes and parameters based on component type
switch compType
    case 'M' % MOSFET: M1 (d g s b) model_name w=... l=... ...
        % Extract nodes (should be in parentheses for Spectre)
        nodes = extract_nodes_from_spectre(line);
        if length(nodes) < 3
            warning('MOSFET %s: insufficient nodes specified (found %d, need at least 3)', name, length(nodes));
            return;
        end
        
        % Ensure we have 4 nodes for MOSFET (drain, gate, source, bulk)
        if length(nodes) == 3
            pins{4} = '0';  % Default bulk to ground if not specified
        end
        
        % Default MOSFET small-signal model parameters
        mosParams = struct('gm', sym('gm'), 'gmb', 0, 'ro', sym('ro'), ...
            'cgs', sym('cgs'), 'cgd', 0, 'cgb', 0, ...
            'cdb', 0, 'csb', 0);
        
        sca = sca.add_mos(sca, name, nodes, mosParams);
        
    case 'R' % Resistor: R1 (n1 n2) resistor r=1k
        nodes = extract_nodes_from_spectre(line);
        if length(nodes) < 2
            warning('Resistor %s: insufficient nodes specified', name);
            return;
        end
        
        % Extract resistance value (look for r=value)
        Rval = extract_parameter_value(line, 'r');
        if isempty(Rval)
            Rval = sym('R'); % Default symbolic value
        end
        
        sca = sca.add_res(sca, name, nodes(1:2), Rval);
        
    case 'C' % Capacitor: C1 (n1 n2) capacitor c=1p
        nodes = extract_nodes_from_spectre(line);
        if length(nodes) < 2
            warning('Capacitor %s: insufficient nodes specified', name);
            return;
        end
        
        % Extract capacitance value (look for c=value)
        Cval = extract_parameter_value(line, 'c');
        if isempty(Cval)
            Cval = sym('C'); % Default symbolic value
        end
        
        sca = sca.add_cap(sca, name, nodes(1:2), Cval);
        
    case 'V' % Voltage Source
        if ignore_sources
            fprintf('Ignoring voltage source: %s\n', name);
            return;
        end
        
        nodes = extract_nodes_from_spectre(line);
        if length(nodes) < 2
            warning('Voltage source %s: insufficient nodes specified', name);
            return;
        end
        
        Vval = 0;  % Always use 0V as specified
        sca = sca.add_vsrc(sca, name, nodes(1:2), Vval);
        
        
    case 'I' % Current Source
        if ignore_sources
            fprintf('Ignoring current source: %s\n', name);
            return;
        end
        
        nodes = extract_nodes_from_spectre(line);
        if length(nodes) < 2
            warning('Current source %s: insufficient nodes specified', name);
            return;
        end
        
        Ival = extract_parameter_value(line, 'dc');
        if isempty(Ival)
            Ival = sym('I');
        end
        sca = sca.add_isrc(sca, name, nodes(1:2), Ival);
        
    otherwise
        return;
end
end

function nodes = extract_nodes_from_spectre(line, debug_flag)
% Extract node names from Spectre format: component (node1 node2 ...)
% Handles both parentheses and space-separated formats

if nargin < 2
    debug_flag = false;
end

nodes = {};

% Look for parentheses first (standard Spectre format)
parenMatch = regexp(line, '\(\s*([^)]+)\s*\)', 'tokens');
if ~isempty(parenMatch)
    nodeStr = strtrim(parenMatch{1}{1});
    nodes = strsplit(nodeStr);
else
    % Fallback: assume nodes are the tokens after component name
    tokens = strsplit(strtrim(line));
    if length(tokens) > 1
        % For lines like: M1 node1 node2 node3 node4 model w=... l=...
        % Skip component name and look for nodes before parameters
        for i = 2:length(tokens)
            token = tokens{i};
            % Stop at parameters (contains =) or known model names
            if contains(token, '=') || any(strcmpi(token, ...
                    {'nch', 'pch', 'nmos', 'pmos', 'resistor', 'capacitor', 'vsource', 'isource'}))
                break;
            end
            % Add as node if it looks like a node name
            if isValidNodeName(token)
                nodes{end+1} = token;
            end
        end
    end
end

% Clean up node names (remove quotes if present)
for i = 1:length(nodes)
    nodes{i} = strrep(nodes{i}, '"', '');
    nodes{i} = strrep(nodes{i}, '''', '');
    nodes{i} = strtrim(nodes{i});
end

% Debug output
if debug_flag
    fprintf('Line: %s\n', line);
    fprintf('Extracted nodes: [%s]\n', strjoin(nodes, ', '));
end
end

function isValid = isValidNodeName(token)
% Check if token is a valid node name
% Valid node names: alphanumeric, underscore, numbers (including 0)
isValid = ~isempty(token) && ischar(token) && ~isempty(regexp(token, '^[a-zA-Z0-9_]+$', 'once'));
end

function value = extract_parameter_value(line, paramName)
% Extract parameter value from Spectre format: param=value

value = [];

% Create regex pattern for parameter extraction
pattern = sprintf('%s\\s*=\\s*([\\w\\.\\-\\+eE]+[a-zA-Z]*)', paramName);
match = regexp(line, pattern, 'tokens', 'ignorecase');

if ~isempty(match)
    valueStr = match{1}{1};
    
    % Try to convert to number first
    numVal = str2double(valueStr);
    if ~isnan(numVal)
        value = numVal;
    else
        % Keep as symbolic if not a simple number
        value = sym(valueStr);
    end
end
end




function sca = update_component_param(sca, compName, paramName, newValue)

found = false;

for j = 1:numel(sca.components)
    comp = sca.components{j};   % extract struct from cell
    
    if strcmp(comp.name, compName)
        if isfield(comp.params, paramName)
            oldVal = comp.params.(paramName);
            fprintf('Updating %s.%s from %s to %s\n', ...
                compName, paramName, string(oldVal), newValue);
            
            comp.params.(paramName) = str2sym(newValue);  % update
            sca.components{j} = comp;                % put it back
            found = true;
            break;
        else
            error('Parameter "%s" not found in %s.params.', ...
                paramName, compName);
        end
    end
end

if ~found
    error('Component %s not found in sca.components.', compName);
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%





function [H_simplified, info] = typicalValueSimplify(H, typical_values, level)
% TYPICALVALUESIMPLIFY - Simplify symbolic transfer functions using typical values
%
% Syntax: [H_simplified, info] = typicalValueSimplify(H, typical_values, level)
%
% Inputs:
%   H - Symbolic transfer function H(s)
%   typical_values - Structure with typical parameter values, e.g.,
%                    struct('gm_ro', 50, 'RL', 1e4, 'CL', 1e-9)
%   level - Negligibility level factor (default: 10)
%
% Outputs:
%   H_simplified - Simplified transfer function
%   info - Structure containing simplification details


if nargin < 3
    level = 10; % Default negligibility level
end

% Initialize info structure
info = struct();
info.original_H = H;
info.level = level;
info.iterations = 0;
info.terms_removed = {};
info.common_factors_removed = {};

fprintf('Starting Typical Value Symbolic Simplification...\n');
fprintf('Negligibility level: %g\n', level);

% Step 1: Decompose H(s) into numerator and denominator
fprintf('\nStep 1: Decomposing H(s)...\n');
[N, D] = numden(H);
fprintf('Numerator N(s) = %s\n', char(N));
fprintf('Denominator D(s) = %s\n', char(D));

% Main iterative loop
max_iterations = 20;
for main_iter = 1:max_iterations
    fprintf('\n=== Main Iteration %d ===\n', main_iter);
    info.iterations = main_iter;
    
    % Store previous state for convergence check
    N_prev = N;
    D_prev = D;
    
    % Step 2-4: Expand, group, and simplify numerator
    fprintf('\nSimplifying numerator...\n');
    [N_grouped, N_terms] = expandAndGroupDimensional(N, typical_values);
    [N, N_removed] = simplifyPolynomialTypical(N_grouped, typical_values, level);
    
    % Step 2-4: Expand, group, and simplify denominator
    fprintf('\nSimplifying denominator...\n');
    [D_grouped, D_terms] = expandAndGroupDimensional(D, typical_values);
    [D, D_removed] = simplifyPolynomialTypical(D_grouped, typical_values, level);
    
    % Store removed terms info
    if main_iter == 1
        info.terms_removed.numerator = N_removed;
        info.terms_removed.denominator = D_removed;
    else
        info.terms_removed.numerator = [info.terms_removed.numerator, N_removed];
        info.terms_removed.denominator = [info.terms_removed.denominator, D_removed];
    end
    
    % Step 5: Remove common factors
    fprintf('\nRemoving common factors...\n');
    [N, D, common_factors] = removeCommonFactors(N, D);
    if ~isempty(common_factors)
        fprintf('Removed common factors: %s\n', char(common_factors));
        info.common_factors_removed{end+1} = common_factors;
    end
    
    % Check for convergence
    if isequal(N, N_prev) && isequal(D, D_prev)
        fprintf('\nConverged after %d iterations\n', main_iter);
        break;
    end
    
    if main_iter == max_iterations
        fprintf('\nReached maximum iterations (%d)\n', max_iterations);
    end
end

% Step 6: Reconstruct simplified transfer function
fprintf('\nStep 6: Reconstructing simplified H(s)...\n');
H_simplified = N / D;
H_simplified = simplify(H_simplified);

fprintf('\nSimplification complete!\n');
fprintf('Original H = %s\n', char(H));
fprintf('Simplified H = %s\n', char(H_simplified));

% Calculate complexity reduction
original_length = length(char(H));
simplified_length = length(char(H_simplified));
reduction_percent = (1 - simplified_length/original_length) * 100;

info.complexity_reduction = reduction_percent;
fprintf('Complexity reduction: %.1f%%\n', reduction_percent);
end

function [grouped_terms, all_terms] = expandAndGroupDimensional(poly, typical_values)
% Expand polynomial and group terms by powers of s with dimensional analysis

syms s

% Expand the polynomial
expanded = expand(poly);

% Convert to polynomial coefficients
coeff_vec = coeffs(expanded, s);
powers = 0:(length(coeff_vec)-1);

% Initialize grouped terms structure
grouped_terms = struct();
all_terms = {};

% Group terms by powers of s
for i = 1:length(coeff_vec)
    power = powers(i);
    coeff = coeff_vec(i);
    
    % Skip zero coefficients
    if coeff == 0
        continue;
    end
    
    % Decompose coefficient into individual terms if it's a sum
    terms = extractTerms(coeff);
    
    % Apply dimensional grouping and parameter unification
    terms = unifyParameters(terms, typical_values);
    
    power_key = sprintf('s%d', power);
    if ~isfield(grouped_terms, power_key)
        grouped_terms.(power_key) = {};
    end
    
    % Add each term to the appropriate power group
    for j = 1:length(terms)
        grouped_terms.(power_key){end+1} = terms{j};
        all_terms{end+1} = terms{j} * s^power;
    end
end

fprintf('Grouped terms by powers of s:\n');
fields = fieldnames(grouped_terms);
for i = 1:length(fields)
    fprintf('  %s: %d terms\n', fields{i}, length(grouped_terms.(fields{i})));
end
end

function unified_terms = unifyParameters(terms, typical_values)
% Unify parameters of the same type (e.g., treat all gm*ro as same entity)

unified_terms = terms; % Start with original terms

for i = 1:length(terms)
    term = terms{i};
    term_str = char(term);
    
    % Replace all gm*ro combinations with unified symbol
    % This is a simplified approach - in practice, you'd want more sophisticated
    % pattern matching for different parameter combinations
    
    % Example unifications (you can extend this):
    % All gm*ro products -> gm_ro_unified
    if contains(term_str, 'gm') && contains(term_str, 'ro')
        % This is a simplified replacement - you may need more sophisticated parsing
        unified_terms{i} = replacePairs(term, 'gm', 'ro', 'gm_ro_unified');
        
    end
    
    % All capacitors -> C_unified
    if contains(term_str, 'C')
        unified_terms{i} = replaceCapacitors(unified_terms{i});
    end
    
    % All resistors -> R_unified
    if contains(term_str, 'R') && ~contains(term_str, 'ro')
        unified_terms{i} = replaceResistors(unified_terms{i});
    end
end
end

function new_term = replacePairs(term, param1, param2, unified_name)
% Replace parameter pairs with unified symbol
% This function finds all combinations of param1*param2 and replaces them

% Convert to string for pattern matching
term_str = char(term);

% Create patterns to match
patterns = {
    [param1, '*', param2],...       % gm*ro
    [param2, '*', param1],...       % ro*gm
    [param1, param2],...            % gmro (without *)
    [param2, param1]                % rogm (without *)
    };

% Replace all patterns with unified name
new_term_str = term_str;
for i = 1:length(patterns)
    pattern = patterns{i};
    new_term_str = strrep(new_term_str, pattern, unified_name);
end

% Convert back to symbolic
new_term = str2sym(new_term_str);

% Alternative approach using symbolic manipulation
try
    % Get all variables in the term
    vars = symvar(term);
    
    % Find param1 and param2 variables
    param1_vars = vars(contains(string(vars), param1));
    param2_vars = vars(contains(string(vars), param2));
    
    if ~isempty(param1_vars) && ~isempty(param2_vars)
        % Create substitution rules
        subs_old = [];
        subs_new = [];
        
        % Handle all combinations of param1*param2
        for i = 1:length(param1_vars)
            for j = 1:length(param2_vars)
                pair_product = param1_vars(i) * param2_vars(j);
                if has(term, pair_product)
                    subs_old = [subs_old, pair_product];
                    subs_new = [subs_new, sym(unified_name)];
                end
            end
        end
        
        % Apply substitutions
        if ~isempty(subs_old)
            new_term = subs(term, subs_old, subs_new);
        else
            new_term = term;
        end
    else
        new_term = term;
    end
catch
    % If symbolic approach fails, use string replacement result
    new_term = sym(new_term_str);
end
end

function new_term = replaceCapacitors(term)
% Replace all capacitor instances with unified symbol

% Convert to string for analysis
term_str = char(term);

% Find all capacitor variables (assuming they start with 'C')
vars = symvar(term);
cap_vars = [];

for i = 1:length(vars)
    var_str = char(vars(i));
    if startsWith(var_str, 'C') && length(var_str) > 1
        % Check if it's likely a capacitor (C followed by numbers/letters)
        if isstrprop(var_str(2), 'alphanum')
            cap_vars = [cap_vars, vars(i)];
        end
    end
end

% Replace all capacitor variables with unified symbol
if ~isempty(cap_vars)
    subs_old = cap_vars;
    subs_new = repmat(sym('C_unified'), size(cap_vars));
    new_term = subs(term, subs_old, subs_new);
else
    new_term = term;
end

% Additional string-based replacement for common patterns
new_term_str = char(new_term);
common_cap_patterns = {'Cgs', 'Cgd', 'Cdb', 'Csb'};

for i = 1:length(common_cap_patterns)
    if contains(new_term_str, common_cap_patterns{i})
        new_term_str = strrep(new_term_str, common_cap_patterns{i}, 'C_unified');
    end
end

% Convert back to symbolic if changes were made
if ~strcmp(char(new_term), new_term_str)
    try
        new_term = sym(new_term_str);
    catch
        % Keep original if conversion fails
    end
end
end

function new_term = replaceResistors(term)
% Replace all resistor instances with unified symbol

% Convert to string for analysis
term_str = char(term);

% Find all resistor variables (assuming they start with 'R' but not 'ro')
vars = symvar(term);
res_vars = [];

for i = 1:length(vars)
    var_str = char(vars(i));
    if startsWith(var_str, 'R') && ~strcmp(var_str, 'ro') && ...
            ~startsWith(var_str, 'ro_') && length(var_str) > 1
        % Check if it's likely a resistor (R followed by numbers/letters)
        if isstrprop(var_str(2), 'alphanum')
            res_vars = [res_vars, vars(i)];
        end
    end
end

% Replace all resistor variables with unified symbol
if ~isempty(res_vars)
    subs_old = res_vars;
    subs_new = repmat(sym('R_unified'), size(res_vars));
    new_term = subs(term, subs_old, subs_new);
else
    new_term = term;
end

% Additional string-based replacement for common patterns
new_term_str = char(new_term);
common_res_patterns = {'Rload', 'Rin', 'Rout', 'Rbias'};

for i = 1:length(common_res_patterns)
    if contains(new_term_str, common_res_patterns{i})
        new_term_str = strrep(new_term_str, common_res_patterns{i}, 'R_unified');
    end
end

% Convert back to symbolic if changes were made
if ~strcmp(char(new_term), new_term_str)
    try
        new_term = sym(new_term_str);
    catch
        % Keep original if conversion fails
    end
end
end

% Additional helper function for more sophisticated gm*ro replacement
function unified_terms = replaceGmRoProducts(terms)
% Specifically handle gm*ro products with better pattern matching

unified_terms = cell(size(terms));

for i = 1:length(terms)
    term = terms{i};
    
    % Find all gm and ro variables
    vars = symvar(term);
    gm_vars = [];
    ro_vars = [];
    
    for j = 1:length(vars)
        var_str = char(vars(j));
        if startsWith(var_str, 'gm')
            gm_vars = [gm_vars, vars(j)];
        elseif startsWith(var_str, 'ro')
            ro_vars = [ro_vars, vars(j)];
        end
    end
    
    % Create all possible gm*ro combinations and replace them
    new_term = term;
    for gm_idx = 1:length(gm_vars)
        for ro_idx = 1:length(ro_vars)
            gm_ro_product = gm_vars(gm_idx) * ro_vars(ro_idx);
            if has(new_term, gm_ro_product)
                new_term = subs(new_term, gm_ro_product, sym('gm_ro_unified'));
            end
        end
    end
    
    unified_terms{i} = new_term;
end
end

function terms = extractTerms(expr)
% Extract individual terms from a symbolic expression

expr_str = char(expr);

if contains(expr_str, '+') || (contains(expr_str, '-') && ~startsWith(strtrim(expr_str), '-'))
    if has(expr, symvar(expr))
        children_expr = children(expr);
        if ~isempty(children_expr)
            terms = num2cell(children_expr);
        else
            terms = {expr};
        end
    else
        terms = {expr};
    end
else
    terms = {expr};
end

% Convert to symbolic
for i = 1:length(terms)
    if ~isa(terms{i}, 'sym')
        terms{i} = sym(terms{i});
    end
end
end

function [simplified_poly, removed_terms] = simplifyPolynomialTypical(grouped_terms, typical_values, level)
% Simplify a polynomial by removing negligible terms using typical values

syms s
removed_terms = {};
iteration = 0;

while true
    iteration = iteration + 1;
    fprintf('  Iteration %d...\n', iteration);
    
    terms_removed_this_iter = false;
    fields = fieldnames(grouped_terms);
    
    % Process each power group
    for i = 1:length(fields)
        power_key = fields{i};
        terms = grouped_terms.(power_key);
        
        if length(terms) <= 1
            continue; % Skip if only one term or empty
        end
        
        fprintf('    Processing %s group with %d terms...\n', power_key, length(terms));
        
        % Evaluate terms using typical values
        term_magnitudes = evaluateTermsTypical(terms, typical_values);
        [dominant_idx, dominant_mag] = findDominantTermTypical(term_magnitudes);
        
        fprintf('      Dominant term: %s (magnitude: %g)\n', ...
            char(terms{dominant_idx}), dominant_mag);
        
        % Prune negligible terms
        terms_to_remove = [];
        for j = 1:length(terms)
            if j ~= dominant_idx
                if term_magnitudes(j) * level < dominant_mag
                    terms_to_remove = [terms_to_remove, j];
                    removed_terms{end+1} = terms{j};
                    terms_removed_this_iter = true;
                    fprintf('        Removing negligible term: %s (mag: %g)\n', ...
                        char(terms{j}), term_magnitudes(j));
                end
            end
        end
        
        % Remove identified terms
        terms(terms_to_remove) = [];
        grouped_terms.(power_key) = terms;
    end
    
    % Check if we should continue iterating
    if ~terms_removed_this_iter
        fprintf('  No more terms to remove. Stopping.\n');
        break;
    end
end

% Reconstruct polynomial from remaining terms
simplified_poly = reconstructPolynomial(grouped_terms);
end

function magnitudes = evaluateTermsTypical(terms, typical_values)
% Evaluate term magnitudes using typical values

magnitudes = zeros(1, length(terms));
param_names = fieldnames(typical_values);


for i = 1:length(terms)
    term = terms{i};
    
    % Get all parameter substitutions
    subs_vars = [];
    subs_vals = [];
    
    
    % Add parameter substitutions
    for k = 1:length(param_names)
        param_name = param_names{k};
        
        
        % Handle normalized parameters
        if strcmp(param_name, 'gm_ro')
            % Substitute all gm*ro combinations with this value
            subs_vars = [subs_vars, sym('gm_ro_unified')];
            subs_vals = [subs_vals, typical_values.(param_name)];
        else
            subs_vars = [subs_vars, sym(param_name)];
            subs_vals = [subs_vals, typical_values.(param_name)];
        end
    end
    
    try
        % Evaluate term magnitude
        val = double(subs(term, subs_vars, subs_vals));
        magnitudes(i) = abs(val);
    catch
        magnitudes(i) = 1; % Default value if evaluation fails
        fprintf('      Warning: Could not evaluate term %s, using default magnitude\n', char(term));
    end
end
end

function [dominant_idx, dominant_mag] = findDominantTermTypical(magnitudes)
% Find the dominant term (largest magnitude)

[dominant_mag, dominant_idx] = max(magnitudes);
end

function [N_new, D_new, common_factors] = removeCommonFactors(N, D)
% Remove common factors between numerator and denominator

% Find GCD of numerator and denominator
common_factors = gcd(N, D);

if common_factors == 1 || common_factors == sym(1)
    % No common factors
    N_new = N;
    D_new = D;
    common_factors = [];
else
    % Remove common factors
    N_new = simplify(N / common_factors);
    D_new = simplify(D / common_factors);
end
end

function poly = reconstructPolynomial(grouped_terms)
% Reconstruct polynomial from grouped terms

syms s
poly = 0;

fields = fieldnames(grouped_terms);
for i = 1:length(fields)
    power_key = fields{i};
    power = str2double(power_key(2:end));
    terms = grouped_terms.(power_key);
    
    % Sum all terms for this power
    coeff_sum = 0;
    for j = 1:length(terms)
        coeff_sum = coeff_sum + terms{j};
    end
    
    if coeff_sum ~= 0
        poly = poly + coeff_sum * s^power;
    end
end

if poly == 0
    poly = sym(0);
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function expr_out = replace_parallel(expr)
    expr_out = mapSymType(expr, 'arith', @detect_replace_parallel);
end


function out = detect_replace_parallel(subexpr)

    out = subexpr;   % default

    [num, den] = numden(simplify(subexpr));
    num_factors = factor(num);

    % We need at least 2 factors in numerator
    if numel(num_factors) < 2
        return
    end

    % Try every pair/triple/... subset of factors as candidate "parallel set"
    for k = 2:numel(num_factors)
        combs = nchoosek(1:numel(num_factors), k);

        for ci = 1:size(combs,1)
            vars = unique(num_factors(combs(ci,:)));

            % Test the identity: does subexpr / (product of leftover multipliers)
            % equal parallel(vars)?
            leftover = prod(num_factors(setdiff(1:numel(num_factors), combs(ci,:))));

            test_expr = simplify(subexpr / leftover);

            lhs = simplify(1/test_expr);
            rhs = simplify(sum(1./vars));

            if isequal(lhs,rhs)
                args = strjoin(arrayfun(@char,vars,'UniformOutput',false), ',');
                out = leftover * sym(['parallel(', args, ')']);
                return
            end
        end
    end
end











