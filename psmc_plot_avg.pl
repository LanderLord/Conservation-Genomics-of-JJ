#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

# =============================================================================
# psmc_plot_avg.pl  --  merged single-file version
#   step 1: parse .psmc output, rescale time (years) and Ne (x10^4)
#   step 2: rebuild every trajectory as a step function, evaluate on the union
#           of all breakpoints, and take the equally weighted arithmetic mean
#   step 3: plot individuals as thin light step lines + mean as a thick red line
#
#   merged from: psmc_plot.pl v0.2.0 (after lh3, with coordinate export),
#                caculate.pl (cross-sample mean), psmc_plot_draw.pl (re-plot)
# =============================================================================

my $version = "1.0.0";

my %opts = (
    u => 2.5e-8,        # absolute mutation rate per site per generation
    s => 100,           # skip used in data preparation
    g => 25,            # years per generation
    n => 20,            # take the n-th EM iteration (0 = keep best-RI round)
    m => 5,             # minimum iteration allowed in RI selection
    x => 10000,         # x-axis lower limit (years)
    X => 0,             # x-axis upper limit, 0 for auto
    Y => 0,             # y-axis upper limit, 0 for auto
    M => '',            # multiline mode: sample names, separated by , or ;
    N => 0.0,           # false negative rate
    f => "Helvetica,16",# font
    w => 1,             # line width of individual curves
    W => 3,             # line width of the mean curve
    P => "right top",   # key position
    T => '',            # figure title
);

getopts('u:s:g:n:m:LM:x:X:Y:N:f:w:W:P:T:pRG', \%opts);

die <<"END" if (@ARGV < 2);
Usage: perl $0 [options] <out.prefix> <in.psmc> [in2.psmc ...]

Options: -u FLOAT   absolute mutation rate per nucleotide [$opts{u}]
         -s INT     skip used in data preparation [$opts{s}]
         -g INT     number of years per generation [$opts{g}]
         -n INT     take n-th iteration, 0 for best RI [$opts{n}]
         -m INT     minimum number of iteration [$opts{m}]
         -M titles  multiline mode, sample names split by , or ; [null]
         -x FLOAT   minimum years plotted [$opts{x}]
         -X FLOAT   maximum years plotted, 0 for auto [0]
         -Y FLOAT   maximum popsize plotted, 0 for auto [0]
         -N FLOAT   false negative rate [0]
         -f STR     font [$opts{f}]
         -w INT     line width of individual curves [$opts{w}]
         -W INT     line width of the mean curve [$opts{W}]
         -P STR     position of the key [$opts{P}]
         -T STR     figure title [null]
         -L         show the last bin
         -p         convert to PDF (with epstopdf)
         -R         do not remove temporary files
         -G         plot grid
END

my $prefix = shift(@ARGV);          # remaining @ARGV are .psmc files

# ---------------------------------------------------------------------------
# 1. read .psmc, keep the chosen iteration, rescale to years and Ne
# ---------------------------------------------------------------------------
my (@data, $d, $N0, $skip, $Mseg, $Msize, $id, $min_ri, $do_store, $gof, $round);
my (@FN, @nscale, @tscale, @alpha, $dt, @xshift, $N_scale, $fn0);
$skip = 1;                          # guard: no MM skip line in some files
$fn0  = $opts{N};

if ($opts{M}) {                     # modifiers: =f *n :t @a >s
    my @t = split(/[,;]+/, $opts{M});
    for my $x (@t) {
        push(@FN,     ($x =~ /=([^\s=*:@>]+)/) ? $1 : 0.0);
        push(@nscale, ($x =~ /\*([^\s=*:@>]+)/) ? $1 : 1.0);
        push(@tscale, ($x =~ /:([^\s=*:@>]+)/) ? $1 : 0.0);
        push(@alpha,  ($x =~ /\@([^\s=*:@>]+)/) ? $1 : 1.0);
        push(@xshift, ($x =~ />([^\s=*:@>]+)/) ? $1 : 0.0);
    }
    for (0 .. $#alpha) {
        $alpha[$_] = 2 * (2 + $alpha[$_]) / 3.0 / (1 + $alpha[$_]);
    }
}

my $max_intv = 10000;
$id = $do_store = 0;
$gof = 'RI';

while (<>) {
    if (/^MM.*skip:(\d+)/) {
        $skip = $1 * $opts{s};
    } elsif (/^MM.*pattern:(\S+?)[,;]?\s*$/) {
        my ($n_lambda, @pat) = parse_pattern($1);
        my $k = $#pat - 1;
        $k-- while ($k >= 0 && $pat[$k] == $pat[$#pat]);
        $max_intv = $k;
    } elsif (/^MM.*is_decoding:/) {
        $d = \%{$data[$id++]};
        $min_ri = 1e30;
    } elsif (/^RD\s(\S+)/) {
        $Mseg = $Msize = $dt = 0;
        $round = $1;
    } elsif (/^(RI|GF)\s(\S+)/ && $1 eq $gof && $2 !~ /nan|inf/) {
        $do_store = 0;
        if ($round >= $opts{m} && $2 < $min_ri) {
            $min_ri = $2;
            $do_store = 1;
        }
        $do_store = ($round == $opts{n}) ? 1 : 0 if ($opts{n} > 0);
    } elsif ($do_store && /^TR\s(\S+)\s(\S+)/) {
        ($d->{T}, $d->{R}) = ($1 / $skip, $2 / $skip);
        $N_scale = 1.0;
        $N_scale /= 1.0 - $FN[$id - 1]  if (defined $FN[$id - 1]);
        $N_scale /= 1.0 - $fn0;
        $N_scale /= $alpha[$id - 1]     if (defined $alpha[$id - 1]);
        $N0 = $1 / $skip / (4 * $opts{u}) * $N_scale;
        $d->{N0} = $N0;
        $d->{RI} = $min_ri;
    } elsif ($do_store && /^DT\s(\S+)/) {
        $dt = $1;
    } elsif ($do_store && /^RS\s(\d+)\s(\S+)\s(\S+)\s*(.*)$/) {
        # bin_index  t_k  lambda_k  [extra columns are ignored]
        my ($bin, $tk, $lk, $rest) = ($1, $2, $3, $4);
        next if (!defined($opts{L}) && $bin > $max_intv);
        my $s = (defined $nscale[$id - 1]) ? $nscale[$id - 1] : 1.0;
        my $t = (defined $tscale[$id - 1]) ? $tscale[$id - 1] : 0.0;
        my $xshift_v = (defined $xshift[$id - 1]) ? $xshift[$id - 1] : 0.0;
        my @extra = split(/\s+/, $rest // '');
        @{$d->{D}[$bin]} = (
            2 * $N0 * ($tk + $dt) * (1.0 - $t) * $opts{g} + $xshift_v,   # x: years
            $lk * $N0 * $s / 10000,                                      # y: Ne x10^4
            @extra);
        $Mseg = $1 if ($Mseg < $bin);
        $Msize = 2 * $N0 * $tk if ($Msize < 2 * $N0 * $tk);
    } elsif ($do_store && /^PA\s(.*)/) {
        $d->{PAR} = $1;
    } elsif ($do_store && /^\/\//) {
        $d->{Mseg} = $Mseg;
        $d->{Msize} = $Msize;
    }
}

die "No PSMC trajectory was read, check -n / -m or the input file\n" unless (@data);

# ---------------------------------------------------------------------------
# 2. build one table per sample, then the cross-sample step-wise mean
# ---------------------------------------------------------------------------
my @titles = $opts{M} ? split(/[,;]+/, $opts{M}) : ('popsize');
my (@curves, %all_x);

for my $i (0 .. $#data) {
    next unless ($data[$i]{D});
    my (@x, @y);
    for my $q (@{$data[$i]{D}}) {
        next unless (defined $q);
        push(@x, $q->[0] + 0.0);
        push(@y, $q->[1] + 0.0);
        $all_x{$q->[0] + 0.0} = 1;
    }
    next unless (@x);
    my $name = defined($titles[$i]) ? $titles[$i] : "sample_$i";
    push(@curves, {name => $name, x => \@x, y => \@y});
}

die "No trajectory left after parsing\n" unless (@curves);

# common time grid: union of all breakpoints, forced to start at 0
my @grid = sort {$a <=> $b} keys %all_x;
unshift(@grid, 0.0) unless (@grid && $grid[0] == 0.0);

my (@avg_x, @avg_y);
for my $xv (@grid) {
    my ($sum, $cnt) = (0.0, 0);
    for my $c (@curves) {
        # step definition: value of the last breakpoint <= interval start
        my $yv = 0.0;
        for (my $k = $#{$c->{x}}; $k >= 0; --$k) {
            if ($c->{x}[$k] <= $xv) { $yv = $c->{y}[$k]; last }
        }
        $sum += $yv;
        $cnt++;
    }
    push(@avg_x, $xv);
    push(@avg_y, $cnt ? $sum / $cnt : 0.0);
}

# ---------------------------------------------------------------------------
# 3. write coordinates and the plotting table
# ---------------------------------------------------------------------------
open(my $sfh, '>', "${prefix}_all_coords.txt") || die "Cannot write ${prefix}_all_coords.txt: $!";
print $sfh "sample\tx_coord\ty_coord\n";
for my $c (@curves) {
    for my $k (0 .. $#{$c->{x}}) {
        printf $sfh "%s\t%.10f\t%.10f\n", $c->{name}, $c->{x}[$k], $c->{y}[$k];
    }
}
close($sfh);

open(my $afh, '>', "${prefix}_avg_coords.txt") || die "Cannot write ${prefix}_avg_coords.txt: $!";
print $afh "x_coord\taverage_y_coord\n";
for my $k (0 .. $#avg_x) {
    printf $afh "%.10f\t%.10f\n", $avg_x[$k], $avg_y[$k];
}
close($afh);

my (@tmp, @titles_safe);
for my $i (0 .. $#curves) {
    my $fn = sprintf("%s.s%02d.txt", $prefix, $i);
    open(my $fh, '>', $fn) || die "Cannot write $fn: $!";
    for my $k (0 .. $#{$curves[$i]{x}}) {
        printf $fh "%.10f\t%.10f\n", $curves[$i]{x}[$k], $curves[$i]{y}[$k];
    }
    close($fh);
    push(@tmp, $fn);
    (my $s = $curves[$i]{name}) =~ s/["\\]/_/g;
    push(@titles_safe, $s);
}
my $avg_file = "$prefix.mean.txt";
open(my $mfh, '>', $avg_file) || die "Cannot write $avg_file: $!";
for my $k (0 .. $#avg_x) {
    printf $mfh "%.10f\t%.10f\n", $avg_x[$k], $avg_y[$k];
}
close($mfh);

# ---------------------------------------------------------------------------
# 4. gnuplot
# ---------------------------------------------------------------------------
my $yran = ($opts{Y} > 0) ? $opts{Y} : '*';
my $xran = ($opts{X} > 0) ? $opts{X} : '*';
my $grid = $opts{G} ? "set grid" : "unset grid";
my $afont = qq/font "$opts{f}"/;
my $ylab_mu = sprintf("%.1fx10^{-8}", $opts{u} / 1e-8);
my $xlab = "Years (g=$opts{g}, {/Symbol m}=$ylab_mu)";
my $ylab = "Effective population size (x10^4)";

open(my $gp, '>', "$prefix.gp") || die "Cannot write $prefix.gp: $!";
print $gp "set size 1, 0.8\n";
print $gp "set xrange [$opts{x}:$xran]\n";
print $gp "set log x\n";
print $gp qq(set format x "10^{\\%L}"\n);
print $gp "set mxtics 10\n";
print $gp "set mytics 10\n";
print $gp "$grid\n";
print $gp "set key $opts{P}\n";
print $gp "set xtics $afont\n";
print $gp "set ytics nomirror $afont\n";
print $gp "set xlabel \"$xlab\" $afont\n";
print $gp "set ylabel \"$ylab\" $afont\n";
print $gp "set terminal postscript enhanced color eps \"$opts{f}\"\n";
print $gp "set output \"$prefix.eps\"\n";
print $gp "set style line 101 lt 1 lc rgb \"#FFA0A0\" lw $opts{w}\n";
print $gp "set style line 200 lt 1 lc rgb \"#FF0000\" lw $opts{W}\n";
print $gp "set title \"$opts{T}\"\n" if ($opts{T});
print $gp "plot ";
my @plots;
for my $i (0 .. $#tmp) {
    push(@plots, sprintf("\"%s\" u 1:2 t \"%s\" w steps ls 101", $tmp[$i], $titles_safe[$i]));
}
push(@plots, sprintf("\"%s\" u 1:2 t \"Average\" w steps ls 200", $avg_file));
print $gp join(", ", @plots), "\n";
close($gp);

system("gnuplot $prefix.gp") == 0
    or warn "gnuplot failed, the script file $prefix.gp is kept for inspection\n";
system("epstopdf $prefix.eps") if (defined $opts{p});

# ---------------------------------------------------------------------------
# 5. clean up
# ---------------------------------------------------------------------------
unless (defined $opts{R}) {
    unlink(@tmp);
    unlink($avg_file);
    unlink("$prefix.gp");
}

printf STDERR "done: %d samples, %d intervals, mean over %d individuals\n",
    scalar(@curves), scalar(@avg_x), scalar(@curves);

# ---------------------------------------------------------------------------
sub parse_pattern {
    my $p = shift;
    $p =~ s/\s+//g;
    my @seg = split(/\+/, $p);
    my $n_lambda = 0;
    my @stack;
    for my $sg (@seg) {
        my ($x1, $x2) = (1, $sg);
        ($x1, $x2) = ($1, $2) if ($sg =~ /(\d+)\*(\d+)/);
        push(@stack, $x2) for (0 .. $x1 - 1);
        $n_lambda += $x1;
    }
    my @ret;
    for my $i (0 .. $#stack) {
        push(@ret, $i) for (0 .. $stack[$i] - 1);
    }
    return ($n_lambda, @ret);
}
