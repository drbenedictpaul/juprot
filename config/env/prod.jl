using Genie.Configuration

const config = Config(
  server_port = parse(Int, get(ENV, "PORT", "8080")),
  server_host = "0.0.0.0",
  log_level = :info,
  log_to_file = false,
  server_handle_static_files = true,
  path_build = "build",
  path_env = "env",
  path_app = "app",
  path_resources = "resources",
  path_lib = "lib",
  path_helpers = "helpers",
  path_log = "log",
  path_tasks = "tasks",
  path_test = "test",
  path_config = "config"
)

ENV["JULIA_REVISE"] = "off"