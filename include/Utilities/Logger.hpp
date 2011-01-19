/*
 * Copyright (C) 2010 Philipp Wähnert
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __LOGGER_HPP
#define __LOGGER_HPP

#include <ostream>

namespace TensorCalculus {
  
  template<typename CharT, typename Traits>
  struct StreamLogger {
    std::basic_ostream<CharT, Traits>* stream;
    char separator;
    
    StreamLogger(std::basic_ostream<CharT, Traits>& stream, char separator = '\t')
      : stream(&stream), separator(separator) { }
    
    template<typename T1>
    void operator () (T1 x1) const
    {
      *stream << x1 << '\n';
      stream->flush();
    }
    
    template<typename T1, typename T2>
    void operator () (T1 x1, T2 x2) const 
    {
      *stream << x1 << separator
              << x2 << '\n';
      stream->flush();
    }
    
    template<typename T1, typename T2, typename T3>
    void operator () (T1 x1, T2 x2, T3 x3) const 
    {
      *stream << x1 << separator 
              << x2 << separator
              << x3 << '\n';
      stream->flush();
    }

    template<typename T1, typename T2, typename T3, typename T4>
    void operator () (T1 x1, T2 x2, T3 x3, T4 x4) const 
    {
      *stream << x1 << separator 
              << x2 << separator
              << x3 << separator
              << x4 << '\n';
      stream->flush();
    }

  };
  
  template<typename CharT, typename Traits>
  StreamLogger<CharT, Traits> stream_logger(std::basic_ostream<CharT, Traits>& stream, char separator = '\t')
  {
    return StreamLogger<CharT, Traits>(stream, separator);
  }
  
  struct NullLogger {
    template<typename T1>
    void operator () (T1 x1) const { }

    template<typename T1, typename T2>
    void operator () (T1 x1, T2 x2) const { }

    template<typename T1, typename T2, typename T3>
    void operator () (T1 x1, T2 x2, T3 x3) const { }

    template<typename T1, typename T2, typename T3, typename T4>
    void operator () (T1 x1, T2 x2, T3 x3, T4 x4) const { }
  };
  
} // namespace TensorCalculus

#endif // __LOGGER_HPP
