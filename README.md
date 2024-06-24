# PitchEffect~ object for MaxMSP

This repository contains source codes for a MaxMSP External Object that can transform the input audio signal to a new pitch. This object can also be used as a vocoder. The underlying principle is based on the minimum-phase response, similar to the STRAIGHT vocoder, which is the predecessor of the [WORLD vocoder](https://github.com/mmorise/World).

## Features

- **Pitch Shifting:** Transform the pitch of an input audio signal to any specified pitch.
- **Vocoding:** Utilise the object as a vocoder for various signal processing applications.

## Optimization

- The object is optimised for **macOS** and does not work on other operating systems, sorry.

## Disclaimer

The developers take no responsibility for any issues or damages that may occur from using this object.

## License

This project is licensed under the [Unlicense](https://unlicense.org/).

## Installation

- Requires [MaxSDK](https://github.com/Cycling74/max-sdk)
- Clone this repository within the MaxSDK directory and use CMake as with other objects

## Usage

See `example.maxpat` for usage instructions.

Feel free to provide feedback to the project via GitHub.
