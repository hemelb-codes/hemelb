-- This file is part of HemeLB and is Copyright (C)
-- the HemeLB team and/or their institutions, as detailed in the
-- file AUTHORS. This software is provided under the terms of the
-- license in the file LICENSE.
on fileLog(message)
	set cmd to "echo '" & (message as text) & "' >> ~/tmp/hemelb-debug.log"
	do shell script cmd
end fileLog

on run args
	set binary to item 1 of args
	set commandFile to ""
	set startIndex to 2

	if item 2 of args is "-file" then
		set commandFile to item 3 of args
		set startIndex to 4
	end if

	set pIds to {}
	repeat with i from startIndex to length of args
		if item i of args is not "" then
			set pIds to pIds & item i of args
		end if
	end repeat

	set lldbPath to do shell script "which lldb"
	if lldbPath is equal to "" then
		set gdbPath to do shell script "which gdb"
		if gdbPath is equal to "" then
			-- Error - no known debugger!
			return "Cannot find a debugger"
		else
			gdbRun(binary, commandFile, pIds)
		end if
	else
		lldbRun(binary, commandFile, pIds)
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
			else
				tell application "System Events" to tell process "Terminal" to keystroke "t" using command down
				delay 1
				set newTab to do script in front tab of front window
			end if
			copy newTab to the end of tabList
			set custom title of newTab to "Rank " & rank
			set rank to rank + 1
		end repeat
	end tell
	return tabList
end CreateTabs

on lldbRun(binary, extraCommandFile, pIds)
	set debugger to "lldb"
    -- Optional user-supplied commands (e.g. setting breakpoints)
	if extraCommandFile is not "" then
		set debugger to debugger & " -s " & extraCommandFile
	end if

    -- Required start up commands
    set resumeCommandFile to myDir() & "/resume.lldb"
    set debugger to debugger & " -s " & resumeCommandFile

	set tabList to CreateTabs(pIds)
	tell application "Terminal"
		-- run commands in tabs
		repeat with i from 1 to count pIds
			set pId to item i of pIds
			set curTab to item i of tabList
			set cmd to debugger & " -p " & pId
			do script cmd in curTab
		end repeat
		
		tell application "System Events" to tell process "Terminal" to keystroke "}" using command down
	end tell
end lldbRun

on gdbRun(binary, commandFile, pIds)
	set debuggerCommandFile to myDir() & "/resume.gdb"
	set debugger to "gdb -x " & debuggerCommandFile
	if commandFile is not "" then
		set debugger to debugger & " -x " & commandFile
	end if
	
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
