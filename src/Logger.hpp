/*!
 *
 * Boost log
 * Why A seperate header file for logging, to seperate complexity of logging
 *
 */
#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>

namespace scs {
namespace internal {

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;

/*!
 *
 * Currently only console logging is supported
 */
class Logger {
 public:
  Logger() {
    logging::core::get()->set_filter(logging::trivial::severity >=
                                     logging::trivial::trace);
    logging::add_common_attributes();
  }
  src::severity_logger<logging::trivial::severity_level> lg;
};
}
}
#endif  // LOGGER_HPP
