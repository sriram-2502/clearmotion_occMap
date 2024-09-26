function [xlog, rho_log, u_cl] = run_nmpc_with_density(estimated_height_map, bin_params, occ_dist, rbf_weights, centers, sigma, x0, xf)
    % mpc_with_density_control: Solves the MPC problem with density constraints using CasADi
    %
    % Parameters:
    %   weights - RBF weights fitted to the height map
    %   centers - RBF centers for the height map
    %   grid_map - The height map (2D array) for terrain
    %   x0 - Initial position of the vehicle [x, y]
    %   xf - Target position of the vehicle [x, y]
    %   N - Prediction horizon
    %   time_total - Total simulation time
    %   dt - Time step for the simulation
    %   sigma - Spread parameter for the RBFs
    %
    % Returns:
    %   xlog - Log of the vehicle's trajectory (states over time)
    %   rho_log - Log of the density function values at each step
    %   u_cl - Closed-loop control inputs over time

    import casadi.*

    % Set defaults
    set(0, 'DefaultFigureColor', 'w');
    set(0, 'DefaultAxesColor', 'w');
    set(0, 'DefaultLineLineWidth', 2);
    set(0, 'DefaultAxesLineWidth', 1);
    
    % setup colors for plots
    colors = colororder;
    blue = colors(1,:);
    red = colors(2,:);
    yellow = colors(3,:);
    green = colors(5,:);
    obsColor = [.7 .7 .7]; % Obstacle color -> Grey

    %% Setup Parameters
    show_plot = true;
    new_rows = bin_params.new_rows;
    new_cols = bin_params.new_cols;

    time_total = 30; % total time for running MPC simulation
    dt = 0.1; % MPC timestep
    N = 10; % length of MPC horizon

    % States for integrator
    x = SX.sym('x');
    y = SX.sym('y');
    states = [x; y];
    n_states = length(states);

    % Control inputs in density space
    u1 = SX.sym('u1');
    u2 = SX.sym('u2');
    controls = [u1; u2];
    n_controls = length(controls);

    % MPC setup
    Q = 10 * diag([1, 1]);
    R = 1 * diag([1, 1]);
    P_terminal = 10 * diag([1, 1]);
    dt_sim = dt;
    C_t = 0.1;

    xmin = [-inf; -inf];
    xmax = -xmin;
    umin = [-1; -1];
    umax = -umin;
    x_ini = x0;

    % Density function setup
    b_x = query_rbf_height(states', rbf_weights, centers, sigma);
    rho_trav = 1 - b_x;
    rho_trav = Function('rho', {states}, {rho_trav});

    %% Dynamics Setup
    [dx_dt, f, g] = integrator_dynamics(states, controls);
    F = Function('F', {states, controls}, {dx_dt});

    %% CasADi MPC Setup
    % State and control decision variables
    X = SX.sym('X', n_states, N + 1);
    U = SX.sym('U', n_controls, N);
    C = SX.sym('C', N);

    % Parameters for initial state and reference state
    P = SX.sym('P', n_states + N * n_states);

    obj = 0;
    constraints = [];

    % Constraint: Initial condition
    st = X(:, 1);
    constraints = [constraints; st - P(1:n_states)];

    % Compute cost and constraints over the horizon
    for k = 1:N
        st = X(:, k);
        con = U(:, k);

        % Cost function (for tracking or stabilization)
        q_x = (st - P(n_states * k + 1:n_states * k + n_states))' * Q * (st - P(n_states * k + 1:n_states * k + n_states));
        obj = obj + q_x + (con' * R * con);

        % Dynamics constraint: x(k+1) = F(x(k), u(k))
        st_next = X(:, k + 1);
        f_value = F(st, con);
        st_next_euler = st + dt * f_value;
        constraints = [constraints; st_next - st_next_euler];
    end

    % Terminal cost
    st = X(:, N + 1);
    terminal_cost = (st - P(n_states + 1:2 * n_states))' * P_terminal * (st - P(n_states + 1:2 * n_states));
    obj = obj + terminal_cost;

    % Density constraint over the horizon
    for k = 1:N
        st = X(:, k);
        rho = rho_trav(st');
        slack = dt * C(k) * rho;
        constraints = [constraints; slack]; % Simplified density constraint
    end

    % Setup optimization problem
    OPT_variables = [reshape(X, n_states * (N + 1), 1); reshape(U, n_controls * N, 1); reshape(C, N, 1)];
    nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', constraints, 'p', P);

    opts = struct;
    opts.ipopt.max_iter = 100;
    opts.ipopt.print_level = 0;
    opts.print_time = 0;
    opts.ipopt.acceptable_tol = 1e-8;
    opts.ipopt.acceptable_obj_change_tol = 1e-6;

    solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

    % Define constraints bounds
    args = struct;
    %------------------- equality constraints for dyanmics -------------------------------
    args.lbg(1:n_states*(N+1)) = 0; 
    args.ubg(1:n_states*(N+1)) = 0;
    
    %------------------ inequality constraints -------------------------------
    % bounds for CONSTRAINT: divergence constriant from MPC-CDF paper
    args.lbg(n_states*(N+1)+1 : n_states*(N+1)+N) = 0; 
    args.ubg(n_states*(N+1)+1 : n_states*(N+1)+N) = inf; 
    
    % bounds for CONSTRAINT: b(x(k))*rho(k) <= gamma
    % one constraint for entire horizon (sum over horizon)
    % args.lbg(n_states*(N+1)+N+1) = 0; 
    % args.ubg(n_states*(N+1)+N+1) = gamma; 
    
    % ------------------- bounds ----------------------------------------------
    args.lbx(1:n_states:n_states*(N+1),1) = xmin(1); %state x lower bound
    args.ubx(1:n_states:n_states*(N+1),1) = xmax(1); %state x upper bound
    args.lbx(2:n_states:n_states*(N+1),1) = xmin(2); %state y lower bound
    args.ubx(2:n_states:n_states*(N+1),1) = xmax(2); %state y upper bound
    
    args.lbx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = umin(1); %u1 lower bound
    args.ubx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = umax(1); %u1 upper bound
    args.lbx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = umin(2); %u2 lower bound
    args.ubx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = umax(2); %u2 upper bound
    
    args.lbx(n_states*(N+1)+n_controls*N+1:n_states*(N+1)+n_controls*N+N,1) = 0; %C lower bound
    args.ubx(n_states*(N+1)+n_controls*N+1:n_states*(N+1)+n_controls*N+N,1) = inf; %C upper bound

    %% Simulate MPC Controller
    t0 = 0;
    u0 = zeros(N, n_controls);
    X0 = repmat(x0, 1, N + 1)';
    C_0 = repmat(C_t, 1, N);
    tstart = 0;
    tend = dt_sim;

    % Initialize logs
    xlog = [];
    u_cl = [];
    rho_log = [];
    tlog = [];
    xlog(:, 1) = x0;
    rho_log = [rho_log, rho_trav(x0')];
    tlog(1) = t0;

    mpciter = 1;
    w_bar = waitbar(0,'1','Name','Simulating MPC-CDF...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    while norm((x0 - xf), 2) > 1e-2 && mpciter < time_total / dt_sim
        max_iter = time_total/dt_sim;
        current_time = mpciter * dt_sim;
        waitbar(mpciter/max_iter,w_bar,sprintf(string(mpciter)+'/'+string(max_iter)))

        % Set parameters (x0, target, etc.)
        t_predict = current_time + N * dt_sim;
        x_target = [2 * t_predict; 13];  % Example target update logic
        x_ref = generate_reference(x0, x_target, N, dt);
        args.p(1:n_states) = x0;
        args.p(n_states + 1:n_states + N * n_states) = x_ref(:);

        % Initial values for optimization variables
        args.x0 = [reshape(X0', n_states * (N + 1), 1); reshape(u0', n_controls * N, 1); reshape(C_0', N, 1)];
        sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, 'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);

        % Extract control input and state trajectory
        u = reshape(full(sol.x(n_states * (N + 1) + 1:n_states * (N + 1) + n_controls * N))', n_controls, N)';
        C_0 = reshape(full(sol.x(end - N + 1:end))', 1, N);

        % Apply the control
        u_star = u(1, :);
        [t, X] = ode45(@(t, x) full(F(x, u_star)), [tstart, tend], x0);
        x0 = X(end, :)';

        % Logging
        xlog(:, mpciter + 1) = x0;
        rho_log = [rho_log, rho_trav(x0')];
        u_cl = [u_cl; u(1, :)];
        tlog(mpciter + 1) = t0;

        % Shift trajectory for the next iteration
        X0 = reshape(full(sol.x(1:n_states * (N + 1)))', n_states, N + 1)';
        X0 = [X0(2:end, :); X0(end, :)];

        % Update time variables
        tstart = tend;
        tend = tend + dt_sim;
        t0 = tend;
        mpciter = mpciter + 1;
    end
    waitbar_graphics = findall(0,'type','figure','tag','TMWWaitbar');
    delete(waitbar_graphics);

    %% show plots
    if(show_plot)
        figure(111)
        subplot(1,2,1)
        % plot heightmap
        [x, y] = meshgrid(1:new_cols, 1:new_rows);  % Original grid
        x = x + occ_dist;
        surf(x, y, estimated_height_map, 'FaceAlpha', 0.5);
        xlabel('X');
        ylabel('Y');
        zlabel('Height');
        view(2)
        hold on
        
        % plot x-y-z trajecotry
        traj = plot(xlog(1,:), xlog(2,:),'-','LineWidth', 2,'Color',red);
        xlabel('x(m)','interpreter','latex','FontSize',20);
        ylabel('y(m)','interpreter','latex','FontSize',20);
        hold on
        
        % plot start and target 
        plot(x_ini(1), x_ini(2), 'o', 'MarkerSize',10, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;
        if(~tracking)
            plot(xf(1), xf(2), 'o', 'MarkerSize',10, 'MarkerFaceColor',green,'MarkerEdgeColor',green); hold on;
        end
        %setup plots
        axes1 = gca;
        box(axes1,'on');
        axis(axes1,'square');
        
        hold(axes1,'off');
        xlabel('position, $x$ (m)','interpreter','latex', 'FontSize', 20);
        ylabel('position, $y$ (m)','interpreter','latex', 'FontSize', 20);
        xlim([1,28]); ylim([1,28])
        
        %%%%%%%%%%%%%%
        subplot(1,2,2)
        % plot heightmap reconstructed using rbfs
        % Generate a fine grid of coordinates
        [x_fine, y_fine] = meshgrid(linspace(1, cols, 10*cols), linspace(1, rows, 10*rows));
        coordinates_fine = [x_fine(:), y_fine(:)];
        
        % Reconstruct the height map using the RBFs and weights on the fine grid
        num_centers = size(centers, 1);
        rbf_matrix_fine = zeros(size(coordinates_fine, 1), num_centers);
        for i = 1:num_centers
            diff = coordinates_fine - centers(i, :);
            rbf_matrix_fine(:, i) = exp(-sum(diff.^2, 2) / (2 * sigma^2));
        end
        fitted_height_map_fine = reshape(rbf_matrix_fine * rbf_weights, size(x_fine));
        
        surf(x_fine, y_fine, fitted_height_map_fine, 'FaceAlpha', 0.5, 'EdgeColor','none');
        xlabel('X');
        ylabel('Y');
        zlabel('Height');
        view(2)
        hold on
        
        % plot x-y-z trajecotry
        traj = plot(xlog(1,:), xlog(2,:),'-','LineWidth', 2,'Color',red);
        xlabel('x(m)','interpreter','latex','FontSize',20);
        ylabel('y(m)','interpreter','latex','FontSize',20);
        hold on
        
        % plot start and target 
        plot(x_ini(1), x_ini(2), 'o', 'MarkerSize',10, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;
        if(~tracking)
            plot(xf(1), xf(2), 'o', 'MarkerSize',10, 'MarkerFaceColor',green,'MarkerEdgeColor',green); hold on;
        end
        
        %setup plots
        axes1 = gca;
        box(axes1,'on');
        axis(axes1,'square');
        xlim([1,28]); ylim([1,28])
        hold(axes1,'off');
        xlabel('position, $x$ (m)','interpreter','latex', 'FontSize', 20);
        ylabel('position, $y$ (m)','interpreter','latex', 'FontSize', 20);
    end
end
