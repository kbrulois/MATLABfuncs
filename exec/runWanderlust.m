function G = runWanderlust(k, l, s, graph, num_landmarks)

data = csvread(getenv('path2wlInput'), 1);

G.Opts.metric = getenv('wl_distMetric');
G.Opts.k = k;
G.Opts.l = l;
G.Opts.num_graphs = graph;
G.Opts.s = s;
G.Opts.num_landmarks = num_landmarks;
G.Opts.verbose = true;
G.Opts.branch = false;
G.Opts.partial_order = [];
G.Opts.deblur = false;
G.Opts.snn = 0;
G.Opts.ann = false;
G.Opts.voting_scheme = getenv('wl_voting_scheme'); % default
G.Opts.band_sample = true;
G.Opts.flock_landmarks = 2;
G.Opts.search_connected_components = true;
G.Opts.plot_landmark_paths = false;
G.Opts.plot_data = [data(:,1) data(:,2)];
G.Opts.lnn = [];
G.Opts.landmarks = [];
G.Opts.disallow = [];
G.Opts.cell_clusters = [];
G.Opts.end_clusters = [];
G.Opts.plot_debug_branch = false;
G.Opts.kEigs = 4;


G = wanderlust(data, G.Opts)

csvwrite(getenv('path2wlOutput'), G.T)