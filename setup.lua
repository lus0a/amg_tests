ug_load_script("ug_util.lua")
--ug_load_script("solver_util/setup_rsamg.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("util/mac_util.lua")
ug_load_script("util/easy_util.lua")


EnableLUA2VM(true)
SetDebugLevel(debugID.LUACompiler, 4)
amgsetup = {}

ug_load_script("setup_amg.lua")
ug_load_script("setup_disc.lua")
ug_load_script("setup_solve_check.lua")

 
------------------------------------------------------------------------------------------------------------------------------------
-- solver check (random solve)

------------------------------------------------------------------------------------------------------------------------------------
-- misc

function getHome()
	local p=os.getenv("HOME")
	if p == nil then
		return "."
	else
		return p
	end
end

function CreateOutputDirectory(name, mbsize)
	local thePath
	name = "/"..name.."/"
	if GetOperatingSystem() == "apple" then
		if mbsize == nil then mbsize = 2048 end
		thePath = util.createMacRAMDisk("ug4RAMDisk", mbsize)..name		
	else
		thePath = getHome()..name
	end
	CreateDirectory(thePath) 
	return thePath		
end

function LogAllOutput(thePath)
	a = GetLogAssistant()
	a:set_output_process(-1)
	if paramOutputProcess ~= -1 and ProcRank() ~= paramOutputProcess then a:enable_terminal_output(false) end
	a:enable_file_output(true, thePath.."log_p"..ProcRank()..".txt")	
end

