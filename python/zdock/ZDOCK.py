#!/usr/bin/env python3

# Copyright (c) 2019, Arjan van der Velde, Weng Lab
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import re


class Structure(object):

    def __init__(self):
        self._filename = ""
        self._translation = (0.0, 0.0, 0.0)
        self._rotation = (0.0, 0.0, 0.0)

    @property
    def filename(self):
        return self._filename

    @property
    def translation(self):
        return self._translation

    @property
    def rotation(self):
        return self._rotation

    @filename.setter
    def filename(self, f):
        self._filename = f

    @translation.setter
    def translation(self, t):
        if len(t) > 2:
            self._translation = (float(t[0]), float(t[1]), float(t[2]))
        else:
            self._translation = (float(t[0]), float(t[1]), 0.0)

    @rotation.setter
    def rotation(self, r):
        if len(r) > 2:
            self._rotation = (float(r[0]), float(r[1]), float(r[2]))
        else:
            self._rotation = (float(r[0]), float(r[1]), 0.0)

    def __repr__(self):
        return "%s\t%.3f\t%.3f\t%.3f" % \
               (self.filename, self.translation[0], self.translation[1],
                self.translation[2])


class Prediction(object):

    def __init__(self):
        self._filename = ""
        self._translation = (0, 0, 0)
        self._rotation = (0.0, 0.0, 0.0)
        self._score = 0.0
        self._ismzdock = False

    @property
    def translation(self):
        return self._translation

    @property
    def rotation(self):
        return self._rotation

    @property
    def ismzdock(self):
        return self._ismzdock

    @property
    def score(self):
        return self._score

    @translation.setter
    def translation(self, t):
        if len(t) > 2:
            self._translation = (int(t[0]), int(t[1]), int(t[2]))
        else:
            self._translation = (int(t[0]), int(t[1]), 0)

    @rotation.setter
    def rotation(self, r):
        if len(r) > 2:
            self._rotation = (float(r[0]), float(r[1]), float(r[2]))
        else:
            self._rotation = (float(r[0]), float(r[1]), 0.0)

    @ismzdock.setter
    def ismzdock(self, m):
        self._ismzdock = (m and m)  # make it a bool

    @score.setter
    def score(self, s):
        self._score = float(s)

    def __repr__(self):
        if self.ismzdock:
            return "%.6f\t%.6f\t%d\t%d\t%.2f" % \
                   (self.rotation[0], self.rotation[1], self.translation[0],
                    self.translation[1], self.score)
        else:
            return "%.6f\t%.6f\t%.6f\t%d\t%d\t%d\t%.3f" % \
                   (self.rotation[0], self.rotation[1], self.rotation[2],
                    self.translation[0], self.translation[1],
                    self.translation[2], self.score)


class ZDOCK(object):

    def __init__(self, filename, n=None):
        self._receptor = Structure()
        self._ligand = Structure()
        self._predictions = []
        self._boxsize = 0
        self._spacing = 0.0
        self._isswitched = False
        self._ismzdock = False
        self._isfixed = False
        self._version = 0
        self._symmetry = 0
        self._filename = os.path.expanduser(os.path.normpath(filename))
        self._read(n)

    @property
    def receptor(self):
        return self._receptor

    @property
    def ligand(self):
        if self.ismzdock:
            raise Exception("ligand() not supported for M-ZDOCK")
        return self._ligand

    @property
    def npredictions(self):
        return len(self._predictions)

    @property
    def predictions(self):
        return self._predictions

    @property
    def ismzdock(self):
        return self._ismzdock

    @property
    def isfixed(self):
        return self._isfixed

    @property
    def symmetry(self):
        if not self.ismzdock:
            raise Exception("symmetry() not supported for M-ZDOCK")
        return self._symmetry

    @property
    def boxsize(self):
        return self._boxsize

    @property
    def filename(self):
        return self._filename

    @property
    def version(self):
        return self._version

    @property
    def spacing(self):
        return self._spacing

    @property
    def lines(self):
        return self._lines()

    @property
    def isswitched(self):
        return self._isswitched

    @ismzdock.setter
    def ismzdock(self, m):
        self._ismzdock = m and m  # make it bool

    @staticmethod
    def _procline(s):
        x = re.sub("^([^#]*).*$", "\\1", s).rstrip()
        if x:
            return x.split('\t')
        return []

    def _read(self, n=None):
        header = []
        with open(self._filename, 'r') as infile:
            self._predictions = []
            headerdone = False
            linenum = 1
            for line in infile:
                p = Prediction()
                p.ismzdock = False
                line = self._procline(line)
                if not line:
                    continue
                if 7 == len(line):
                    if self.ismzdock:
                        # M-ZDOCK established but 7-column prediction encountered
                        raise Exception("Invalid M-ZDOCK prediction (line " + str(linenum) + ")")
                    # ZDOCK prediction
                    (a, b, c, d, e, f, g) = \
                        [t(s) for t, s in zip((float, float, float, int, int, int, float), line)]
                    p.rotation = (a, b, c)
                    p.translation = (d, e, f)
                    p.score = g
                    self._predictions.append(p)
                    headerdone = True
                elif 5 == len(line):
                    if not self.ismzdock and self.npredictions > 0:
                        # ZDOCK established but 5-column prediction encountered
                        raise Exception("Invalid ZDOCK prediction (line " + str(linenum) + ")")
                    # M-ZDOCK prediction
                    (a, b, c, d, e) = [t(s) for t, s in zip((float, float, int, int, float), line)]
                    p.rotation = (a, b, 0.0)
                    p.translation = (c, d, 0)
                    p.score = e
                    p.ismzdock = True
                    self.ismzdock = True
                    self._predictions.append(p)
                    headerdone = True
                else:
                    if not headerdone:
                        header.append(line)
                    else:
                        raise Exception("Invalid prediction (line " + str(linenum) + ")")
                linenum += 1
                if n is not None and len(self._predictions) > n:
                    break

        # figure out header
        if not self.ismzdock:
            # ZDOCK
            if 5 == len(header):
                self._version = 1  # new style has 5 header rows
                try:
                    # string '0' to bool; str -> int -> bool
                    (self._boxsize, self._spacing, self._isswitched) = \
                        [t(s) for t, s in zip((int, float, lambda x: bool(int(x))), header[0])]
                    self._isfixed = False
                except Exception as e:
                    raise Exception("ZDOCK header error: " + str(e))
            elif 4 == len(header):
                self._version = 0  # old style has 4 header rows
                self._isfixed = True
                self._isswitched = False
                try:
                    (self._boxsize, self._spacing) = [t(s) for t, s in zip((int, float), header[0])]
                except Exception as e:
                    raise Exception("ZDOCK header error: " + str(e))
            else:
                raise Exception("ZDOCK header error; ZDOCK header must have 4 or 5 rows")
        else:
            # M-ZDOCK
            if 3 != len(header):
                raise Exception("M-ZDOCK header error; M-ZDOCK header must have 3 rows")
            try:
                (self._boxsize, self._spacing, self._symmetry) = \
                    [t(s) for t, s in zip((int, float, int), header[0])]
            except Exception as e:
                raise Exception("M-ZDOCK header error: " + str(e))
            if 3 > self.symmetry:
                raise Exception("M-ZDOCK symmetry cannot be less than 3")

        # receptor
        if not self.isfixed:
            try:
                if 3 == len(header[1]):
                    (a, b, c) = [t(s) for t, s in zip((float, float, float), header[1])]
                    self._receptor.rotation = (a, b, c)
                elif 2 == len(header[1]) and self.ismzdock:
                    (a, b) = [t(s) for t, s in zip((float, float), header[1])]
                    self._receptor.rotation = (a, b, 0.0)
                else:
                    raise Exception("Unable to obtain receptor initial rotation")
            except Exception as e:
                raise Exception("Unable to obtain receptor initial rotation: " + str(e))
        try:
            if self.isswitched:
                h = 4 - (not self.version)
            else:
                h = 3 - (not self.version)
            (a, b, c, d) = [t(s) for t, s in zip((str, float, float, float), header[h])]
            self._receptor.filename = a
            self._receptor.translation = (b, c, d)
        except Exception as e:
            raise Exception("Unable to obtain receptor initial translation: " + str(e))

        # ligand
        if not self.ismzdock:
            try:
                if self.isfixed:
                    h = 1
                else:
                    h = 2
                if 3 == len(header[h]):
                    (a, b, c) = [t(s) for t, s in zip((float, float, float), header[h])]
                    self._ligand.rotation = (a, b, c)
                else:
                    raise Exception("Unable to obtain ligand initial rotation")
            except Exception as e:
                raise Exception("Unable to obtain ligand initial rotation: " + str(e))
            try:
                if self.isswitched:
                    h = 3 - (not self.version)
                else:
                    h = 4 - (not self.version)
                (a, b, c, d) = [t(s) for t, s in zip((str, float, float, float), header[h])]
                self._ligand.filename = a
                self._ligand.translation = (b, c, d)
            except Exception as e:
                raise Exception("Unable to obtain ligand initial translation: " + str(e))

    def _lines(self):
        if self.ismzdock:
            yield "%d\t%.1f\t%d" % (self.boxsize, self.spacing, self.symmetry)
            yield "%.6f\t%.6f\t%.6f" % \
                  (self._receptor.rotation[0], self._receptor.rotation[1], self._receptor.rotation[2])
        elif self.isfixed:
            yield "%d\t%.1f" % (self.boxsize, self.spacing)
            yield "%.6f\t%.6f\t%.6f" % \
                  (self._ligand.rotation[0], self._ligand.rotation[1], self._ligand.rotation[2])
        else:
            yield "%d\t%.1f\t%d" % (self.boxsize, self.spacing, self.isswitched)
            yield "%.6f\t%.6f\t%.6f" % self._receptor.rotation
            yield "%.6f\t%.6f\t%.6f" % self._ligand.rotation
        if not self.ismzdock:
            if self.isswitched:
                yield str(self._ligand)
                yield str(self._receptor)
            else:
                yield str(self._receptor)
                yield str(self._ligand)
        else:
            yield str(self._receptor)
        for p in self._predictions:
            yield str(p)

    def __repr__(self):
        return '\n'.join([x for x in self._lines()])
