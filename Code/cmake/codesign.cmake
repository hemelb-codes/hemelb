# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

# On macOS you need to sign your executable to avoid popups and to
# allow a debugger to attach.
#
# 1. Create a certificate through EITHER XCode OR KeyChain.
#   I. XCode
#     a) Open XCode
#     b) Go to preferences -> accounts and make sure your Apple ID is
#        there.
#     c) Click "Manage Certificates", if you don't have one, click the
#        "+" and add an "Apple Development" certificate.
#
#   II. KeyChain
#     a) Open KeyChain Access
#     b) Keychain Access > Certificate Assistant > Create a certificate
#     c) Choose: Name: something without spaces for easy command line
#        use. ID type: self signed root. Cert. type: code signing.
#     d) Make sure you enable "Let me override defaults" and continue.
#     e) Serial and validity as you like (note name + serial must be
#        unique).
#     f) Your name etc can be whatever.
#     g) Key info: default are fine. Keep hitting continue until you get
#        to the "keychain", then choose system.
#
# 2. In KeyChain Access, trust this certificate for code signing.
#
# 3. (Optionally) Trust codesign to access the certificate's private
#    key without a password. Expand the triangle by the certifcate's
#    name to show the key. Double click and go to "access control",
#    choose the "+" to add an app, and select /usr/bin/codesign.
#
# 4. Tell CMake to use the identity by setting the CODESIGN_IDENTITY
#    variable.
#

if (APPLE)
  set(_default_cs ON)
  set(_default_entitlements
    "${CMAKE_CURRENT_LIST_DIR}/macos-allow-debug.plist")
else()
  set(_default_cs OFF)
  set(_default_entitlements "")
endif()

define_property(TARGET PROPERTY CODESIGN INHERITED
  BRIEF_DOCS "Do code-signing"
  FULL_DOCS "Whether to sign the executable")
option(CODESIGN "Default for whether to code-sign executables" ${_default_cs})

define_property(TARGET PROPERTY CODESIGN_IDENTITY INHERITED
  BRIEF_DOCS "Code-signing identity name (macOS only)"
  FULL_DOCS "Name of identity to use for code-signing")
set(CODESIGN_IDENTITY ""
  CACHE STRING
  "Default name of identity to use for code-signing")

define_property(TARGET PROPERTY CODESIGN_ENTITLEMENTS INHERITED
  BRIEF_DOCS "Entitlements file (macOS only)"
  FULL_DOCS "Entitlements plist file to use when code-signing executables")
set(CODESIGN_ENTITLEMENTS "${_default_entitlements}"
  CACHE STRING
  "Default entitlements file to use for code-signing")

# Helper that first checks the requested property before falling back
# to the global var.
macro(_get_cs_property varname target propname)
  get_target_property("${varname}" "${target}" "${propname}")
  if (NOT "${${varname}}" OR "${${varname}}" STREQUAL "${varname}-NOTFOUND")
    # fall back to global var
    set("${varname}" "${${propname}}")
  endif()
endmacro()

# Warn user if the specified codesigning ID is not available.
#
# Argument is the name of a variable to fill with TRUE/FALSE if we
# found the identity or not.
function(_codesign_check_identity identity_name output_varname)
  if ("$ENV{_codesign_found_${identity_name}}")
    set(${output_varname} TRUE PARENT_SCOPE)
    return()
  endif()

  # Get list of valid codesigning identities from security
  execute_process(
    COMMAND security find-identity -v -p codesigning
    RESULT_VARIABLE rc
    OUTPUT_VARIABLE lines
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  if (rc EQUAL 0)
    # Look for the footer including preceeding newline and whitespace.
    set(nIdRe "\n *([0-9]+) valid identities found$")
    if (lines MATCHES "${nIdRe}")
      # Get number of valid IDs
      set(NVALID ${CMAKE_MATCH_1})
      # Remove footer
      string(REGEX REPLACE "${nIdRe}" "" lines "${lines}")
      # Convert to ;-list
      string(REPLACE "\n" ";" lines "${lines}")
      # Sanity check
      list(LENGTH lines len)
      if (NOT len EQUAL NVALID)
	message(WARNING "security output parse failure")
	set(${output_varname} FALSE PARENT_SCOPE)
	return()
      endif()

      # Split each line `<i>) <hexdigits> "<identity name>"`
      foreach(line "${lines}")
	separate_arguments(items UNIX_COMMAND "${line}")
	# Sanity check
	list(LENGTH items nel)
	if (NOT nel EQUAL 3)
	  message(WARNING "security output parse failure")
	  set(${output_varname} FALSE PARENT_SCOPE)
	  return()
	endif()
	list(GET items 2 id_name)
	if (id_name STREQUAL identity_name)
	  # We are OK!
	  set(ENV{_codesign_found_${identity_name}} TRUE)
	  message(STATUS "Found codesigning identity '${identity_name}'")
	  set(${output_varname} TRUE PARENT_SCOPE)
	  return()
	endif()
      endforeach()
      message(WARNING "Could not find identity '${identity_name}'")
    else()
      message(WARNING "security output not in expected format")
    endif()
  else()
    message(WARNING "Error trying to check codesigning identity")
  endif()
  set(${output_varname} FALSE PARENT_SCOPE)
endfunction()

function(codesign target)
  _get_cs_property(do_cs ${target} CODESIGN)
  if(do_cs)
    find_program(
      CODESIGN_EXECUTABLE NAMES codesign
      DOC "Codesigning executable"
      )

    _get_cs_property(id "${target}" CODESIGN_IDENTITY)
    if (id)
      _codesign_check_identity("${id}" id_found)
      if (id_found)
	# Force re-signing
	set(cs_args "-f")
	list(APPEND cs_args -s "${id}")
	#set(cmd "${cmd} -s ${id}")
	_get_cs_property(entitlements "${target}" CODESIGN_ENTITLEMENTS)
	if (entitlements)
	  list(APPEND cs_args --entitlements "${entitlements}")
	  #set(cmd "${cmd} --entitlements ${entitlements}")
	endif()

	list(APPEND cs_args "$<TARGET_FILE:${target}>")
	#set(cmd "${cmd} $<TARGET_FILE:${target}>")
	add_custom_command(
          TARGET "${target}" POST_BUILD
          COMMAND "${CODESIGN_EXECUTABLE}"
	  ARGS ${cs_args}
          COMMENT "Code-signing ${target} (may prompt for your password)"
          VERBATIM
	  )
      endif()
    endif()
  endif()
endfunction()
