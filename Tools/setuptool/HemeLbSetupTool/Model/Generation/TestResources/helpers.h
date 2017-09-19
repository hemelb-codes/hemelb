// -*- mode: C++; -*-
#ifndef HEMELBSETUPTOOL_TESTRESOURCES_HELPERS_H
#define HEMELBSETUPTOOL_TESTRESOURCES_HELPERS_H

template<class... Types>
const std::string GetResource(const Types&... paths);

template<>
const std::string GetResource() {
  if (const char* env_dir = std::getenv("TESTRESOURCES")) {
    return std::string(env_dir);
  } else {
    throw std::runtime_error("Must set TESTRESOURCES env var");
  }
}

template<class T>
const std::string GetResource(const T& pth) {
  return GetResource() + "/" + pth;
}

template<class T, class... Types>
const std::string GetResource(const T& pth, const Types&... paths) {
  return GetResource(pth) + "/" + GetResource(paths...);
}
#endif
