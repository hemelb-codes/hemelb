on fileLog(message)
	set cmd to "echo '" & (message as text) & "' >> ~/tmp/hemelb-debug.log"
	do shell script cmd
end fileLog

on run args
	set binary to item 1 of args
	
	set pIds to {}
	repeat with i from 2 to length of args
		if item i of args is not "" then
			set pIds to pIds & item i of args
		end if
	end repeat
	
	set lldbPath to do shell script "which lldb"
	if lldbPath is equal to ""
		set gdbPath to do shell script "which gdb"
		if gdbPath is equal to ""
			-- Error - no known debugger!
			return "Cannot find a debugger"
		else
			gdbRun(binary, pIds)
		end if
	else
		lldbRun(binary, pIds)
	end if
end run

on myDir()
	set savedDelimiters to AppleScript's text item delimiters
	try
		set AppleScript's text item delimiters to "/"
		set answer to (text items 1 thru -2 of the POSIX path of (path to me)) as string
		set AppleScript's text item delimiters to savedDelimiters
	on error m number n
		set AppleScript's text item delimiters to savedDelimiters
		error m number n
	end try
	return answer
end myDir

on CreateTabs(pIds)
	set tabList to {}
	tell application "Terminal"
		activate
		-- create tabs
		set rank to 0
		repeat with pId in pIds
			if rank is equal to 0 then
				set newTab to do script
				set wId to (id of first window whose first tab is newTab)
			else
				tell application "System Events" to tell process "Terminal" to keystroke "t" using command down
				set newTab to do script in window id wId
			end if
			copy newTab to the end of tabList
			set custom title of newTab to "Rank " & rank
			set rank to rank + 1
		end repeat
	end tell
	return tabList
end CreateTabs

on lldbRun(binary, pIds)
	set debugger to "lldb"
		
	set tabList to CreateTabs(pIds)
	tell application "Terminal"
		-- run commands in tabs
		repeat with i from 1 to count pIds
			set pId to item i of pIds
			set curTab to item i of tabList
			set cmd to debugger & " -f " & binary & " -p " & pId
			set newTab to do script cmd in curTab
			delay 0.5
			do script ("frame select -r 3") in curTab
			do script ("expr amWaiting = 0") in curTab
			do script ("breakpoint set -F hemelb::debug::ActiveDebugger::BreakHere()") in curTab
			do script ("breakpoint command add -o 'finish' 1") in curTab
			do script ("continue") in curTab
		end repeat
		
		tell application "System Events" to tell process "Terminal" to keystroke "}" using command down
	end tell
end lldbRun

on gdbRun(binary, pIds)
	set debuggerCommandFile to myDir() & "/resume.gdb"
	set debugger to "gdb -x " & debuggerCommandFile
	
	set tabList to CreateTabs(pIds)
	tell application "Terminal"
		-- run commands in tabs
		repeat with i from 1 to count pIds
			set pId to item i of pIds
			set curTab to item i of tabList
			set cmd to debugger & " " & binary & " " & pId
			set newTab to do script cmd in curTab
		end repeat
		
		tell application "System Events" to tell process "Terminal" to keystroke "}" using command down
	end tell
end gdbRun
