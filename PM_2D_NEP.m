%**************************************************************************
%*****************************PM-ABL v2.0**********************************
%**************************************************************************
%**************************************************************************
% 2-D Langevin model to calculate ensemble of particle trajectories
% from U and V winds, variances, and relaxation timescales
%--------------------------------------------------------------------------
%instead of repeated diurnal cycles, the trajectories are run for 48 hours
%with the injection occuring during the first 24 hours in the current setup

%inversion heights are read in from ZINV for the aerosol concentration
%calculation

%C_0 parameter changes as a function of precipitation. Case is considered
%precipitating if the surface rain rate is greater than 0.01 mm/day
%--------------------------------------------------------------------------
case_name = 0220180808;

% Parameters
ens_num = 3; %number of ships
part_num = 1000; %number of particles
x_range = 500; % km
y_range = 500; % km
%particle decay timescale 
tau_decay = 2; %days
%C_0 for precipitating cases
C_0_prec = 0.15;
%C_0 for non-precipitating cases
C_0_noprec = 0.5;
dt = 900; % time step of particle model in seconds
injection_rate = 10^16; %particles per second
rho_mg = 1.15*10^6; %g/m3
total_deleted_from_first_group = 0.0;
% Define the injection schedule
injection_duration = 24 * 3600; % hours --> seconds
total_cycle_duration = 48 * 3600; % hours --> seconds

% Time parameters
day_s = 86400; %number of seconds in day
num_days = 2; %number of diurnal cycles to be simulated
time_s = double((time - time(1)) * 3600.0 * 24.0);
dt_les = round(time_s(2) - time_s(1));
start_t = floor(min(time_s));
end_t = time_s(end);
%create num_days long array of dt_les
time_les = transpose(0:dt_les:end_t);
new_time = transpose(start_t:dt:end_t);
new_time_hr = new_time / 3600.0;
num_step = 1;

%plotting parameters
save_interval = 900; %in seconds
density_gif_initialized = false;
prob_gif_initialized = false;
save_size = ((day_s * num_days) / save_interval) + 1;

% Interpolate U_bl, U2_bl, V_bl, V2_bl, TLu_bl, and TLv_bl
MEAN_U_interp = interp1(time_les, U_bl, new_time);
MEAN_U2_interp = interp1(time_les, U2_bl, new_time);
MEAN_V_interp = interp1(time_les, V_bl, new_time);
MEAN_V2_interp = interp1(time_les, U2_bl, new_time);
TS_u_interp = interp1(time_les, TLu_bl, new_time).*3600.0; %convert to s
TS_v_interp = interp1(time_les, TLv_bl, new_time).*3600.0; %convert to s

% Interpolate precipitation array
prec_thresh = 0.01; %prec. threshold in mm/day
Prec_interp = interp1(time_les, PREC, new_time);

% Interpolate inversion height for particle density calculation
Bl_depth = interp1(time_les, ZINV, new_time).*1000.0; %height in meters

%Interpolate NAc 
NAc_interp = interp1(time_les, NAc_bl, new_time);

% Initialize arrays to store particle positions and velocities
part_pos_x = [];
part_pos_y = [];
part_vel_u = [];
part_vel_v = [];
part_release_time = [];

% Initialize particle positions based on first LES datapoint (t = 2.0 hr)
std_dev_x = 1.0 * 1000.0; % initial x-plume width
std_dev_y = 1.0 * 1000.0; % initial y-plume width

% Initialize ship positions randomly within the domain
initial_positions_x = x_range * 1000 * rand(ens_num, 1); % random x positions in meters
initial_positions_y = y_range * 1000 * rand(ens_num, 1); % random y positions in meters

% Random ship motions
ship_motion_y_options = [-10, 0, 10];
ship_motion_x_options = [-10, 0, 10];
[ship_motion_x_grid, ship_motion_y_grid] = meshgrid(ship_motion_x_options, ship_motion_y_options);
valid_ship_motion_indices = find(~(ship_motion_x_grid == 0 & ship_motion_y_grid == 0));

% Select random indices for ship motions
num_valid_motions = length(valid_ship_motion_indices);
random_indices = randi(num_valid_motions, ens_num, 1); % Generate random indices, allowing repetition
selected_ship_motion_indices = valid_ship_motion_indices(random_indices);
ship_motion_x = ship_motion_x_grid(selected_ship_motion_indices);
ship_motion_y = ship_motion_y_grid(selected_ship_motion_indices);

% Initialize particle velocities
std_dev_vel = 0.0;
mean_u_vel = MEAN_U_interp(1);
mean_v_vel = MEAN_V_interp(1);

% Define GIF parameters
density_filename = 'particle_density.gif';
prob_filename = 'prob_conc.gif';
age_filename = 'age_particle.gif';
gif_filename = 'particle_positions.gif';

% Define grid parameters for density calculation
grid_size = x_range ; % 1 x 1 km grid
x_edges = linspace(0, x_range * 1000, grid_size + 1);
y_edges = linspace(0, y_range * 1000, grid_size + 1);

x_range_m = x_range * 1000; % Convert to meters
y_range_m = y_range * 1000; % Convert to meters

block_size = 5; % 5-km coarse graining
x_blocks = floor((x_range) / block_size);
y_blocks = floor((y_range) / block_size);

%preallocate memory for saved 2D concentrations
particle_conc_2D = zeros(x_blocks, y_blocks, save_size);

% Main time loop to calculate each ensemble particle trajectory
for n = 1:length(new_time)
    current_time = new_time(n);

    % Determine the current time within the total cycle duration
    time_in_cycle = mod(current_time, total_cycle_duration);

    if time_in_cycle <= injection_duration
        % Inject new particles at the ship's current position for each ensemble member
        for m = 1:ens_num
            % Update ship position based on its motion
            initial_positions_x(m) = initial_positions_x(m) + ship_motion_x(m) * dt;
            initial_positions_y(m) = initial_positions_y(m) + ship_motion_y(m) * dt;

            % Check if the ship has left the domain
            if initial_positions_x(m) < 0 || initial_positions_x(m) > x_range * 1000 || ...
                    initial_positions_y(m) < 0 || initial_positions_y(m) > y_range * 1000
                % Regenerate the ship's position and motion
                initial_positions_x(m) = x_range * 1000 * rand();
                initial_positions_y(m) = y_range * 1000 * rand();

                % Select new random motion for the ship
                random_index = randi(num_valid_motions, 1);
                selected_ship_motion_index = valid_ship_motion_indices(random_index);
                ship_motion_x(m) = ship_motion_x_grid(selected_ship_motion_index);
                ship_motion_y(m) = ship_motion_y_grid(selected_ship_motion_index);
            end

            part_pos_init_x = initial_positions_x(m) + std_dev_x * randn(part_num, 1); % particle x position in meters
            part_pos_init_y = initial_positions_y(m) + std_dev_y * randn(part_num, 1); % particle y position in meters
            part_vel_u_init = mean_u_vel + std_dev_vel * randn(part_num, 1);
            part_vel_v_init = mean_v_vel + std_dev_vel * randn(part_num, 1);

            part_pos_x = [part_pos_x; part_pos_init_x];
            part_pos_y = [part_pos_y; part_pos_init_y];
            part_vel_u = [part_vel_u; part_vel_u_init];
            part_vel_v = [part_vel_v; part_vel_v_init];
            part_release_time = [part_release_time; repmat(current_time / 3600, part_num, 1)]; % release time in hours
        end
    end

    % Update particle positions and velocities for all particles
    for i = 1:length(part_pos_x)

        %check for precipitation to determine C_0 value
        if Prec_interp(n) > prec_thresh
            % Calculate relaxation timescale from k/eps from LES
            TLu = TS_u_interp(n) * (1 / (0.75 * C_0_prec));
            TLv = TS_v_interp(n) * (1 / (0.75 * C_0_prec));
        else
            TLu = TS_u_interp(n) * (1 / (0.75 * C_0_noprec));
            TLv = TS_v_interp(n) * (1 / (0.75 * C_0_noprec));
        end

        % Zonal velocities and positions
        part_vel_u(i) = part_vel_u(i) + (MEAN_U_interp(n) - part_vel_u(i)) * (dt / TLu) ...
            + sqrt((2 * MEAN_U2_interp(n) * dt) / TLu) * randn();

        part_pos_x(i) = part_pos_x(i) + part_vel_u(i) * dt;

        % Apply periodic boundary conditions in x direction
        if part_pos_x(i) < 0
            part_pos_x(i) = part_pos_x(i) + x_range * 1000; % wrap to the right side of the domain
        elseif part_pos_x(i) > x_range * 1000
            part_pos_x(i) = part_pos_x(i) - x_range * 1000; % wrap to the left side of the domain
        end

        % Meridional velocities and positions
        part_vel_v(i) = part_vel_v(i) + (MEAN_V_interp(n) - part_vel_v(i)) * (dt / TLv) ...
            + sqrt((2 * MEAN_V2_interp(n) * dt) / TLv) * randn();

        part_pos_y(i) = part_pos_y(i) + part_vel_v(i) * dt;

        % Apply periodic boundary conditions in y direction
        if part_pos_y(i) < 0
            part_pos_y(i) = part_pos_y(i) + y_range * 1000; % wrap to the top side of the domain
        elseif part_pos_y(i) > y_range * 1000
            part_pos_y(i) = part_pos_y(i) - y_range * 1000; % wrap to the bottom side of the domain
        end

    end

    %Track remaining particles after applying decay
    remaining_particles = true(1, length(part_release_time));
    
    %determine unique release times
    unique_release_times = unique(part_release_time);

    for rt = length(unique_release_times):-1:1 % Loop from the last to the first to avoid indexing issues
        release_time = unique_release_times(rt);
        indices = find(part_release_time == release_time);

        particle_age = current_time - release_time * 3600.0; % age in seconds
        decay_prob = 1 - exp(-particle_age / (tau_decay * day_s)); % decay probability based on age

        % Number of particles remaining after applying decay probability
        remaining_fraction = 1 - decay_prob;

        % Update the remaining particles based on the remaining fraction
        if ~isempty(indices)
            num_particles_remaining = round(remaining_fraction * (ens_num * part_num));
            num_particles_to_delete = length(indices) - num_particles_remaining;

            if num_particles_to_delete > 0
                % Randomly select particles to delete
                particles_to_delete = indices(randperm(length(indices), num_particles_to_delete));
                remaining_particles(particles_to_delete) = false; % Mark for deletion
            end
        end

    end

    % Delete the particles after marking
    part_pos_x(~remaining_particles) = [];
    part_pos_y(~remaining_particles) = [];
    part_vel_u(~remaining_particles) = [];
    part_vel_v(~remaining_particles) = [];
    part_release_time(~remaining_particles) = [];

    % Print remaining particles after deletion
    disp(['Number of particles after deletion: ', num2str(length(part_release_time))]);

    % Plot the positions every save_interval
    if mod(current_time, save_interval) == 0
        
        % Calculate and plot particle density/concentration
        figure('visible', 'off');
        %raw particle count per grid cell
        particle_density = histcounts2(part_pos_x, part_pos_y, x_edges, y_edges);

        dx = x_edges(2) - x_edges(1);
        dy = y_edges(2) - y_edges(1); 
        cell_area = dx * dy;  % area of grid cell in m^2

        count_per_particle = (injection_rate * dt) / part_num; %# particles per superparticle
        part_dens = (particle_density * count_per_particle) / cell_area; %# of particles / m2
        particle_concentration = part_dens ./ (Bl_depth(n) * rho_mg); % # / mg 

         %Coarse grain the particle_concentration matrix
        block_size = 5;
        [x_size, y_size] = size(particle_concentration);
        x_blocks = floor(x_size / block_size);
        y_blocks = floor(y_size / block_size);

        coarse_particle_concentration = zeros(x_blocks, y_blocks);

        for i = 1:x_blocks
            for j = 1:y_blocks
                x_start = (i-1)*block_size + 1;
                x_end = i*block_size;
                y_start = (j-1)*block_size + 1;
                y_end = j*block_size;

                % Extract the block and calculate the mean
                block = particle_concentration(x_start:x_end, y_start:y_end);
                coarse_particle_concentration(i, j) = mean(block(:));
            end
        end

        %save 2D particle concentration
        particle_conc_2D(:,:,num_step) = coarse_particle_concentration;
        num_step = num_step + 1; 

        cmap = colormap('jet'); 
        cmap(1, :) = [0 0 0];

        imagesc('XData', x_edges / 1000, 'YData', y_edges / 1000, 'CData', coarse_particle_concentration');
        set(gca, 'FontSize', 10, 'FontName', 'Georgia');
        colorbar;
        colormap(cmap);
        caxis([0 100]); 
        xlim([0 x_range])
        ylim([0 y_range])
        xlabel('X (km)','FontName', 'Georgia');
        ylabel('Y (km)','FontName', 'Georgia');
        title(['Particle Concentration (#/mg) at Time = ', num2str(current_time / 86400), ' Day']);
        grid on;

        % Capture the density plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        % Write to the density GIF file
        if ~density_gif_initialized
            imwrite(imind, cm, density_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
            density_gif_initialized = true;
        else
            imwrite(imind, cm, density_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
        end

        close(gcf);

        concentration_values = coarse_particle_concentration(:);
        concentration_values = concentration_values(~isnan(concentration_values) & ~isinf(concentration_values));
        
        % Calculate the histogram of concentrations
        bin_size = 2;
        min_conc = 0;
        max_conc = 100;
        bin_edges = min_conc:bin_size:max_conc;
        [counts, ~] = histcounts(concentration_values, bin_edges);
        
        % Calculate the bin centers for plotting
        bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
        total_counts = sum(counts);
        probabilities = counts / total_counts;

        % Plot the PDF
        figure('visible', 'off');
        plot(bin_centers, probabilities, 'LineWidth', 4, 'Color', 'black');
        xlabel('Aerosol Concentration (#/mg)');
        ylabel('Probability');
        title(['Probability Distribution of Aerosol Concentration at Time = ', num2str(current_time / 86400), ' Day']);
        ylim([0 0.3]);
        xlim([0 100]);
        grid on;

        % Capture the plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        % Write to the PDF GIF file
        if ~prob_gif_initialized % Initialize the GIF file
            imwrite(imind, cm, prob_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
            prob_gif_initialized = true;
        else
            imwrite(imind, cm, prob_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end

        close(gcf);

    end

end


%------------------------------------------------------
%write particle model data to matlab file
%------------------------------------------------------

%save .mat files
save(['part_conc_' num2str(case_name) '_' num2str(ens_num) 'ship.mat'],'particle_conc_2D');

