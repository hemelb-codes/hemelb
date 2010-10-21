# Get out of sleep and disable the waiting flag
up
up
up
set var amWaiting = 0

# Set up a breakpoint for our break function and a command to finish
# it automatically
break hemelb::debug::ActiveDebugger::BreakHere
commands
finish
end

# carry on
continue
