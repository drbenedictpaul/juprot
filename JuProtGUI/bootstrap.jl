(pwd() != @__DIR__) && cd(@__DIR__) # allow starting app from bin/ dir

using JuProtGUI
const UserApp = JuProtGUI
JuProtGUI.main()
