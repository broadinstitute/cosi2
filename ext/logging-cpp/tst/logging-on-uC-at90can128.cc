/*******************************************************************************
 *
 * Copyright (c) 2008, 2009 Michael Schulze <mschulze@ivs.cs.uni-magdeburg.de>
 * All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without
 *    modification, are permitted provided that the following conditions
 *    are met:
 *
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *
 *    * Neither the name of the copyright holders nor the names of
 *      contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * $Id: logging-on-uC-at90can128.cc 1 2010-01-11 14:24:28Z perch $
 *
 ******************************************************************************/

#define LOGGING_DEFINE_EXTENDED_OUTPUT_TYPE
#include "logging/logging.h"
using namespace ::logging;

#include <avr/io.h>

template<uint32_t baudrate, uint32_t CPU_FREQ>
struct UARTat90can128 {
    UARTat90can128 () {
        uint16_t value=((CPU_FREQ/8/baudrate)-1)/2;
        UBRR0H = (unsigned char) (value>>8);
        UBRR0L = (unsigned char) value;

        /* Enable UART transmitter */
        UCSR0B = ( 1 << TXEN0 );

        /* Set frame format: 8N1 */
        UCSR0C = (1<<UCSZ01)|(1<<UCSZ00);

    }

    UARTat90can128 & operator<<(const char c) {
         /* Wait for empty transmit buffer */
         while ( !(UCSR0A & (1<<UDRE0)) );
         /* Start transmittion */
         UDR0 = c;
         return *this;
    }
};

typedef ::logging::OutputLevelSwitchDisabled <
            ::logging::OutputStream <
                UARTat90can128<19200,16000000>
            >
        > UartLogType;

LOGGING_DEFINE_OUTPUT( UartLogType )

// logging levels can be disabled at compile time
//LOGGING_DISABLE_LEVEL(::logging::Error);
//LOGGING_DISABLE_LEVEL(::logging::Trace);
//LOGGING_DISABLE_LEVEL(::logging::Warning);
//LOGGING_DISABLE_LEVEL(::logging::Info);

int main(int, char**) {
    ::logging::log::emit() << "Hello World! with the logging framework"
        << ::logging::log::endl << ::logging::log::endl;
    ::logging::log::emit() << "Print 15 in hexadecimal "
        << ::logging::log::hex << 15 << ::logging::log::endl;
    ::logging::log::emit() << "Print 15 in decimal "
        << ::logging::log::dec << 15 << ::logging::log::endl;
    ::logging::log::emit() << "Print 15 in octal "
        << ::logging::log::oct << 15 << ::logging::log::endl;
    ::logging::log::emit() << "Print 15 in binary with a tab"
        << ::logging::log::bin << ::logging::log::tab << 15
        << ::logging::log::endl << ::logging::log::endl;

    ::logging::log::emit< ::logging::Error>() << "Logging an Error"
        << ::logging::log::endl;
    ::logging::log::emit< ::logging::Trace>() << "Logging a Trace"
        << ::logging::log::endl;
    ::logging::log::emit< ::logging::Warning>() << "Logging a Warning"
        << ::logging::log::endl;
    ::logging::log::emit< ::logging::Info>() << "Logging an Info"
        << ::logging::log::endl;
    return 0;
}
