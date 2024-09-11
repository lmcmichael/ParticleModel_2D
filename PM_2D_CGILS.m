% Simplified Langevin model to calculate ensemble of particle trajectories
% from U and V winds, variances, and relaxation timescales

% Parameters
ens_num = 3;
part_num = 1000;
x_range = 500; % km
y_range = 500; % km
%particle decay timescale 
tau_decay = 2; %days
%C_0 for precipitating cases
C_0 = 0.15;
%C_0 for non-precipitating cases
%C_0 = 0.5;
dt = 900; % time step of particle model in seconds
injection_rate = 10^16; %particles per second
bl_depth = 1000.0; %boundary layer depth in meters
rho_mg = 1.15*10^6; %g/m3
total_deleted_from_first_group = 0.0;
% Define the injection schedule
injection_duration = 24 * 3600; % hours --> seconds
total_cycle_duration = 24 * 3600; % hours --> seconds

% Time parameters
day_s = 86400; %number of seconds in day
num_days = 11; %number of diurnal cycles to be simulated
time_s = double((time - 12.0) * 3600.0);
dt_les = round(time_s(2) - time_s(1));
start_t = floor(min(time_s));
end_t = num_days * day_s;
%create num_days long array of dt_les
time_les = transpose(0:dt_les:end_t);
time0 = time - 12.0;
t = linspace(0.5, 36, 72);
t_s(:) = t(:) * 3600.0; % s
time_start_index = 1; % (time_start_index/2.0 = # of hours)
t_start = fix(time_start_index);
dt_data = 900.0; % time step in input data (s)
new_time = transpose(start_t:dt:end_t);
new_time_hr = new_time / 3600.0;
num_step = 1;

%plotting parameters
save_interval = 21600; %in seconds (6 hours: 21600 s)
density_gif_initialized = false;
prob_gif_initialized = false;
age_gif_initialized = false;
save_size = (day_s / save_interval) + 1;

%copy diurnal cycle from day 1 repeatedly 
day_index = find(new_time==day_s);
MEAN_U_day = MEAN_U(1:day_index);
MEAN_U2_day = MEAN_U2(1:day_index);
MEAN_V_day = MEAN_U(1:day_index);
MEAN_V2_day = MEAN_V2(1:day_index);
NEW_TS_day = NEW_TS(1:day_index);

%duplicate array for num_days
MEAN_U_dup = repmat(MEAN_U_day, num_days, 1);
MEAN_U2_dup = repmat(MEAN_U2_day, num_days, 1);
MEAN_V_dup = repmat(MEAN_V_day, num_days, 1);
MEAN_V2_dup = repmat(MEAN_V2_day, num_days, 1);
NEW_TS_dup = repmat(NEW_TS_day, num_days, 1);

%clip to time array
MEAN_U_1 = MEAN_U_dup(1:length(new_time));
MEAN_U2_1 = MEAN_U2_dup(1:length(new_time));
MEAN_V_1 = MEAN_V_dup(1:length(new_time));
MEAN_V2_1 = MEAN_V2_dup(1:length(new_time));
NEW_TS_1 = NEW_TS_dup(1:length(new_time));

% Interpolate MEAN_U, MEAN_U2, MEAN_V, MEAN_V2, and NEW_TS
MEAN_U_interp = interp1(time_les, MEAN_U_1, new_time);
MEAN_U2_interp = interp1(time_les, MEAN_U2_1, new_time);
MEAN_V_interp = interp1(time_les, MEAN_V_1, new_time);
MEAN_V2_interp = interp1(time_les, MEAN_V2_1, new_time);
NEW_TS_interp = interp1(time_les, NEW_TS_1, new_time);

% Initialize arrays to store particle positions and velocities
part_pos_x = [];
part_pos_y = [];
part_vel_u = [];
part_vel_v = [];
part_release_time = [];

% Initialize particle positions based on first LES datapoint (t = 2.0 hr)
std_dev_x = 0.5 * width_1std_100_control(t_start) * 1000.0; % 1 std dev in meters
std_dev_y = 100.0; % small y-standard deviation

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

% Convert from two hours to new time index
dt_ind = fix((double(t_s(t_start)) / dt_data) + 1);

% Initialize particle velocities
std_dev_vel = 0.0;
mean_u_vel = MEAN_U(dt_ind);
mean_v_vel = MEAN_U(dt_ind);

% Calculate first time point in new time array
dt_ind_start = round((double(dt_data) / double(dt)) * (dt_ind - 1));

% Define GIF parameters
density_filename = 'particle_density_random_ships.gif';
prob_filename = 'prob_conc.gif';
age_filename = 'age_particle.gif';

% Define grid parameters for density calculation
grid_size = x_range ; % 1 x 1 km grid
x_edges = linspace(0, x_range * 1000, grid_size + 1);
y_edges = linspace(0, y_range * 1000, grid_size + 1);

x_range_m = x_range * 1000; % Convert to meters
y_range_m = y_range * 1000; % Convert to meters

% Coarse grain the particle_concentration_young matrix
block_size = 5; % 5-km coarse graining
x_blocks = floor((x_range) / block_size);
y_blocks = floor((y_range) / block_size);

%preallocate memory for saved 2D concentrations
particle_conc_2D = zeros(x_blocks, y_blocks, save_size);
particle_conc_2D_young = zeros(x_blocks, y_blocks, save_size);
particle_conc_2D_mid = zeros(x_blocks, y_blocks, save_size);
particle_conc_2D_old = zeros(x_blocks, y_blocks, save_size);

peak_concentration = zeros(1, length(new_time));

% Main time loop to calculate each ensemble particle trajectory
for n = dt_ind_start:length(new_time)
    current_time = new_time(n);

    % Determine the current time within the 24-hour cycle
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
            part_release_time = [part_release_time; repmat(current_time / 3600, part_num, 1)]; % release time in hoursm
        end
    end

    % Update particle positions and velocities for all particles
    for i = 1:length(part_pos_x)
        % Calculate relaxation timescale from k/eps from LES
        T_L = NEW_TS_interp(n) * (1 / (0.75 * C_0));

        % Zonal velocities and positions
        part_vel_u(i) = part_vel_u(i) + (MEAN_U_interp(n) - part_vel_u(i)) * (dt / T_L) ...
            + sqrt((2 * MEAN_U2_interp(n) * dt) / T_L) * randn();

        part_pos_x(i) = part_pos_x(i) + part_vel_u(i) * dt;

        % Apply periodic boundary conditions in x direction
        if part_pos_x(i) < 0
            part_pos_x(i) = part_pos_x(i) + x_range * 1000; % wrap to the right side of the domain
        elseif part_pos_x(i) > x_range * 1000
            part_pos_x(i) = part_pos_x(i) - x_range * 1000; % wrap to the left side of the domain
        end

        % Meridional velocities and positions
        part_vel_v(i) = part_vel_v(i) + (MEAN_V_interp(n) - part_vel_v(i)) * (dt / T_L) ...
            + sqrt((2 * MEAN_V2_interp(n) * dt) / T_L) * randn();

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
    day_save = 10.0; %save the day 10 data
    if mod(current_time, save_interval) == 0

         % Calculate particle ages
        particle_ages = (current_time - part_release_time * 3600.0) / day_s; % Age in days

        % Categorize particles based on age
        young_particles = particle_ages <= 1; % 0-1 days
        middle_particles = particle_ages > 1 & particle_ages <= 3; % 1-3 days
        old_particles = particle_ages > 3; % 3+ days

        % Create a 3-panel figure
        figure('visible', 'off');

        % --- Left Panel: 0-1 Day Particles Concentration ---
        % Calculate the 1-km binned concentration for young particles
        particle_density_young = histcounts2(part_pos_x(young_particles), part_pos_y(young_particles), x_edges, y_edges) / (dt * part_num);
        count_per_particle_young = (injection_rate * dt) / part_num; %# particles per particle
        part_dens_young = particle_density_young .* count_per_particle_young; % # of injected particles / m2
        particle_concentration_young = part_dens_young ./ (bl_depth * rho_mg); % # / mg 

        % Coarse grain the particle_concentration_young matrix
        block_size = 5; % 5-km coarse graining
        [x_size, y_size] = size(particle_concentration_young);
        x_blocks = floor(x_size / block_size);
        y_blocks = floor(y_size / block_size);

        coarse_particle_concentration_young = zeros(x_blocks, y_blocks);

        for i = 1:x_blocks
            for j = 1:y_blocks
                x_start = (i-1)*block_size + 1;
                x_end = i*block_size;
                y_start = (j-1)*block_size + 1;
                y_end = j*block_size;

                % Extract the block and calculate the mean
                block = particle_concentration_young(x_start:x_end, y_start:y_end);
                coarse_particle_concentration_young(i, j) = mean(block(:));
            end
        end

        % Plot the coarse-grained concentration for 0-1 day particles
        subplot(1, 3, 1);
        cmap = colormap('jet'); 
        cmap(1, :) = [0 0 0]; % Set the first color (for zero concentration) to black

        imagesc('XData', x_edges / 1000, 'YData', y_edges / 1000, 'CData', coarse_particle_concentration_young');
        set(gca, 'YDir', 'normal');
        colorbar;
        colormap(cmap);
        caxis([0 200]); 
        axis equal; 
        xlim([0 x_range])
        ylim([0 y_range])
        xlabel('X Position (km)');
        ylabel('Y Position (km)');
        title('0-1 Day Aerosol Concentration','FontSize',8);
        grid on;

        % --- Middle Panel: 1-3 Day Particles Concentration ---
        % Calculate the 1-km binned concentration for middle-aged particles
        particle_density_middle = histcounts2(part_pos_x(middle_particles), part_pos_y(middle_particles), x_edges, y_edges) / (dt * part_num);
        count_per_particle_middle = (injection_rate * dt) / part_num; %# particles per particle
        part_dens_middle = particle_density_middle .* count_per_particle_middle; % # of injected particles / m2
        particle_concentration_middle = part_dens_middle ./ (bl_depth * rho_mg); % # / mg 

        % Coarse grain the particle_concentration_middle matrix
        coarse_particle_concentration_middle = zeros(x_blocks, y_blocks);

        for i = 1:x_blocks
            for j = 1:y_blocks
                x_start = (i-1)*block_size + 1;
                x_end = i*block_size;
                y_start = (j-1)*block_size + 1;
                y_end = j*block_size;

                % Extract the block and calculate the mean
                block = particle_concentration_middle(x_start:x_end, y_start:y_end);
                coarse_particle_concentration_middle(i, j) = mean(block(:));
            end
        end

        % Plot the coarse-grained concentration for 1-3 day particles
        subplot(1, 3, 2);
        imagesc('XData', x_edges / 1000, 'YData', y_edges / 1000, 'CData', coarse_particle_concentration_middle');
        set(gca, 'YDir', 'normal');
        colorbar;
        colormap(cmap);
        caxis([0 200]); 
        axis equal; 
        xlim([0 x_range])
        ylim([0 y_range])
        xlabel('X Position (km)');
        ylabel('Y Position (km)');
        title('1-3 Day Aerosol Concentration','FontSize',8);
        grid on;

        % --- Right Panel: 3+ Day Particles Concentration ---
        % Calculate the 1-km binned concentration for old particles
        particle_density_old = histcounts2(part_pos_x(old_particles), part_pos_y(old_particles), x_edges, y_edges) / (dt * part_num);
        count_per_particle_old = (injection_rate * dt) / part_num; %# particles per particle
        part_dens_old = particle_density_old .* count_per_particle_old; % # of injected particles / m2
        particle_concentration_old = part_dens_old ./ (bl_depth * rho_mg); % # / mg 

        % Coarse grain the particle_concentration_old matrix
        coarse_particle_concentration_old = zeros(x_blocks, y_blocks);

        for i = 1:x_blocks
            for j = 1:y_blocks
                x_start = (i-1)*block_size + 1;
                x_end = i*block_size;
                y_start = (j-1)*block_size + 1;
                y_end = j*block_size;

                % Extract the block and calculate the mean
                block = particle_concentration_old(x_start:x_end, y_start:y_end);
                coarse_particle_concentration_old(i, j) = mean(block(:));
            end
        end

        % Plot the coarse-grained concentration for 3+ day particles
        subplot(1, 3, 3);
        imagesc('XData', x_edges / 1000, 'YData', y_edges / 1000, 'CData', coarse_particle_concentration_old');
        set(gca, 'YDir', 'normal');
        colorbar;
        colormap(cmap);
        caxis([0 200]); 
        axis equal; 
        xlim([0 x_range])
        ylim([0 y_range])
        xlabel('X Position (km)');
        ylabel('Y Position (km)');
        title('3+ Day Aerosol Concentration','FontSize',8);
        grid on;

        % Capture the combined plot as an image
        frame_age = getframe(gcf);
        im_age = frame2im(frame_age);
        [imind_age, cm_age] = rgb2ind(im_age, 256);

        % Write to the 3-panel GIF file
        if ~age_gif_initialized % Initialize the GIF file for age categories
            imwrite(imind_age, cm_age, age_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
            age_gif_initialized = true;
        else
            imwrite(imind_age, cm_age, age_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end

        close(gcf);


        % Calculate and plot particle density/concentration
        figure('visible', 'off');
        particle_density = histcounts2(part_pos_x, part_pos_y, x_edges, y_edges) / (dt * part_num);
        count_per_particle = (injection_rate * dt) / part_num; %# particles per particle
        part_dens = particle_density .* count_per_particle; % # of injected particles / m2
        particle_concentration = part_dens ./ (bl_depth * rho_mg); % # / mg 

         % Coarse grain the particle_concentration matrix
        block_size = 5;
        [x_size, y_size] = size(particle_concentration);
        x_blocks = floor(x_size / block_size);
        y_blocks = floor(y_size / block_size);

        coarse_particle_concentration = zeros(x_blocks, y_blocks);

        % Compute the peak concentration for this time step
        peak_concentration(n) = max(particle_concentration(:));

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

        if current_time >= day_save*day_s && current_time <= (day_save+1.0)*day_s 
            particle_conc_2D_young(:,:,num_step) = coarse_particle_concentration_young;
            particle_conc_2D_mid(:,:,num_step) = coarse_particle_concentration_middle;
            particle_conc_2D_old(:,:,num_step) = coarse_particle_concentration_old;
            particle_conc_2D(:,:,num_step) = coarse_particle_concentration;
        
            num_step = num_step + 1;
        end

        cmap = colormap('jet'); 
        cmap(1, :) = [0 0 0];

        imagesc('XData', x_edges / 1000, 'YData', y_edges / 1000, 'CData', coarse_particle_concentration');
        set(gca, 'YDir', 'normal');
        colorbar;
        colormap(cmap);
        caxis([0 200]); % Set color axis limits
        xlim([0 x_range])
        ylim([0 y_range])
        xlabel('X Position (km)');
        ylabel('Y Position (km)');
        title(['Particle Concentration (#/mg) at Time = ', num2str(current_time / 86400), ' Day']);
        grid on;

        % Capture the density plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        % Write to the density GIF file
        if ~density_gif_initialized
            imwrite(imind, cm, density_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
            density_gif_initialized = true;
        else
            imwrite(imind, cm, density_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
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
 
        if current_time >= day_save*day_s && current_time <= (day_save+1.0)*day_s 
            % Save the probabilities array
            save(['probabilities_', num2str(n), '.mat'], 'probabilities');
            save(['bincenters_', num2str(n), '.mat'], 'bin_centers');
        end

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
save(['particle_conc_young_' num2str(ens_num) 'ship.mat'], 'particle_conc_2D_young');
save(['particle_conc_mid_' num2str(ens_num) 'ship.mat'], 'particle_conc_2D_mid');
save(['particle_conc_old_' num2str(ens_num) 'ship.mat'], 'particle_conc_2D_old');
save(['particle_conc_2D_' num2str(ens_num) 'ship.mat'],'particle_conc_2D');


