%% simulation data init
function struct = init_struct_analyse(w,Nts)
    struct.w = w;
    struct.Nts = Nts;

    struct.Ts = 5.0000e-04;
    struct.Tend = 2;

    struct.forward = 1;

    struct.noise_amp = 1e-4;
    struct.init_error_amp = 1e-3;
    struct.init_param_error_amp = 1e-4;

    struct.model = 'tokamak';

    % struct.weight = 1*ones(1,struct.w);
    struct.weight = 1*linspace(0.1,1,struct.w);

    % struct.dweight = 1e-1*ones(1,struct.w);
    struct.dweight = 1e-2*linspace(1,0.1,struct.w);

    struct.ddweight = 0*ones(1,struct.w);

    struct.scale_factor = [1 1e-4 1];
end