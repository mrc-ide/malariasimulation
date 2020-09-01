/*
 * Trompeloeil C++ mocking framework
 *
 * Copyright Björn Fahller 2014-2019
 * Copyright Tore Martin Hagen 2019
 *
 *  Use, modification and distribution is subject to the
 *  Boost Software License, Version 1.0. (See accompanying
 *  file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 *
 * Project home: https://github.com/rollbear/trompeloeil
 */


#ifndef TROMPELOEIL_TESTTHAT_HPP_
#define TROMPELOEIL_TESTTHAT_HPP_

#include "../trompeloeil.hpp"

#ifndef TESTTHAT_HPP
#error "<testthat.h> must be included before <testthat/trompeloeil.hpp>"
#endif

namespace trompeloeil
{
  template <>
  inline void reporter<specialized>::send(
    severity s,
    const char* file,
    unsigned long line,
    const char* msg)
  {
    std::ostringstream os;
    if (line) os << file << ':' << line << '\n';
    os << msg;
    auto failure = os.str();
    if (s == severity::fatal)
    {
      CATCH_FAIL(failure);
    }
    else
    {
      CATCH_CAPTURE(failure);
      CATCH_CHECK(failure.empty());
    }
  }

  template <>
  inline void reporter<specialized>::sendOk(
    const char* trompeloeil_mock_calls_done_correctly)
  {      
      CATCH_REQUIRE(trompeloeil_mock_calls_done_correctly != 0);
  }
}


#endif //TROMPELOEIL_TESTTHAT_HPP_
