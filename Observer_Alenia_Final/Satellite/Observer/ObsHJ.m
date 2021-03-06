function HJ = ObsHJ(Chi11,Chi12,Chi13,Chi21,Chi22,Chi23,Chi31,Chi32,Chi33,Chi41,Chi42,Chi43,Wgps,Wsigma,Wuwb,d_1_2,d_1_3,d_1_4,d_2_3,d_2_4,d_3_4,zg_1,zg_2,zg_3,zg_4,zu_1_2,zu_1_3,zu_1_4,zu_2_3,zu_2_4,zu_3_4)
%OBSHJ
%    HJ = OBSHJ(CHI11,CHI12,CHI13,CHI21,CHI22,CHI23,CHI31,CHI32,CHI33,CHI41,CHI42,CHI43,WGPS,WSIGMA,WUWB,D_1_2,D_1_3,D_1_4,D_2_3,D_2_4,D_3_4,ZG_1,ZG_2,ZG_3,ZG_4,ZU_1_2,ZU_1_3,ZU_1_4,ZU_2_3,ZU_2_4,ZU_3_4)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    01-Jun-2020 11:02:42

t2 = Chi11-Chi21;
t3 = Chi12-Chi22;
t4 = Chi13-Chi23;
t5 = t2.^2;
t6 = t3.^2;
t7 = t4.^2;
t8 = t5+t6+t7;
t9 = Chi11-Chi31;
t10 = Chi12-Chi32;
t11 = Chi13-Chi33;
t12 = t9.^2;
t13 = t10.^2;
t14 = t11.^2;
t15 = t12+t13+t14;
t16 = Chi11-Chi41;
t17 = Chi12-Chi42;
t18 = Chi13-Chi43;
t19 = t16.^2;
t20 = t17.^2;
t21 = t18.^2;
t22 = t19+t20+t21;
t24 = Chi11.*2.0;
t27 = Chi21.*2.0;
t23 = t24-t27;
t31 = Chi31.*2.0;
t25 = -t24+t31;
t35 = Chi41.*2.0;
t26 = -t24+t35;
t28 = t23.^2;
t29 = sqrt(t8);
t30 = d_1_2-t29;
t32 = t24-t31;
t33 = sqrt(t15);
t34 = d_1_3-t33;
t36 = t24-t35;
t37 = sqrt(t22);
t38 = d_1_4-t37;
t39 = 1.0./sqrt(t8);
t40 = 1.0./t8;
t41 = Wuwb.*t28.*t40.*zu_1_2.*(1.0./2.0);
t42 = 1.0./t8.^(3.0./2.0);
t43 = Wuwb.*t28.*t30.*t42.*zu_1_2.*(1.0./2.0);
t44 = 1.0./sqrt(t15);
t45 = t32.^2;
t46 = 1.0./t15;
t47 = 1.0./t15.^(3.0./2.0);
t48 = Wuwb.*t34.*t45.*t47.*zu_1_3.*(1.0./2.0);
t49 = 1.0./sqrt(t22);
t50 = t36.^2;
t51 = 1.0./t22;
t52 = 1.0./t22.^(3.0./2.0);
t53 = Wuwb.*t38.*t50.*t52.*zu_1_4.*(1.0./2.0);
t54 = Wsigma.*2.0;
t55 = Wgps.*zg_1.*2.0;
t56 = Wuwb.*t30.*t39.*zu_1_2.*2.0;
t57 = Wuwb.*t34.*t44.*zu_1_3.*2.0;
t58 = Wuwb.*t38.*t49.*zu_1_4.*2.0;
t60 = Chi12.*2.0;
t63 = Chi22.*2.0;
t59 = t60-t63;
t65 = Chi32.*2.0;
t61 = -t60+t65;
t67 = Chi42.*2.0;
t62 = -t60+t67;
t64 = t59.^2;
t66 = t60-t65;
t68 = t60-t67;
t69 = Wuwb.*t40.*t64.*zu_1_2.*(1.0./2.0);
t70 = Wuwb.*t30.*t42.*t64.*zu_1_2.*(1.0./2.0);
t71 = t66.^2;
t72 = Wuwb.*t34.*t47.*t71.*zu_1_3.*(1.0./2.0);
t73 = t68.^2;
t74 = Wuwb.*t38.*t52.*t73.*zu_1_4.*(1.0./2.0);
t76 = Chi13.*2.0;
t79 = Chi23.*2.0;
t75 = t76-t79;
t81 = Chi33.*2.0;
t77 = -t76+t81;
t83 = Chi43.*2.0;
t78 = -t76+t83;
t80 = t75.^2;
t82 = t76-t81;
t84 = t76-t83;
t85 = Wuwb.*t40.*t80.*zu_1_2.*(1.0./2.0);
t86 = Wuwb.*t30.*t42.*t80.*zu_1_2.*(1.0./2.0);
t87 = t82.^2;
t88 = Wuwb.*t34.*t47.*t87.*zu_1_3.*(1.0./2.0);
t89 = t84.^2;
t90 = Wuwb.*t38.*t52.*t89.*zu_1_4.*(1.0./2.0);
t91 = -t41-t43+t56;
t92 = Chi21-Chi31;
t93 = Chi22-Chi32;
t94 = Chi23-Chi33;
t95 = t92.^2;
t96 = t93.^2;
t97 = t94.^2;
t98 = t95+t96+t97;
t99 = Chi21-Chi41;
t100 = Chi22-Chi42;
t101 = Chi23-Chi43;
t102 = t99.^2;
t103 = t100.^2;
t104 = t101.^2;
t105 = t102+t103+t104;
t106 = t27-t31;
t107 = t27-t35;
t108 = t106.^2;
t109 = sqrt(t98);
t110 = d_2_3-t109;
t111 = t107.^2;
t112 = sqrt(t105);
t113 = d_2_4-t112;
t114 = 1.0./sqrt(t98);
t115 = 1.0./t98;
t116 = Wuwb.*t108.*t115.*zu_2_3.*(1.0./2.0);
t117 = 1.0./t98.^(3.0./2.0);
t118 = Wuwb.*t108.*t110.*t117.*zu_2_3.*(1.0./2.0);
t119 = 1.0./sqrt(t105);
t120 = 1.0./t105;
t121 = Wuwb.*t111.*t120.*zu_2_4.*(1.0./2.0);
t122 = 1.0./t105.^(3.0./2.0);
t123 = Wuwb.*t111.*t113.*t122.*zu_2_4.*(1.0./2.0);
t124 = t56-t69-t70;
t125 = Wgps.*zg_2.*2.0;
t126 = Wuwb.*t110.*t114.*zu_2_3.*2.0;
t127 = Wuwb.*t113.*t119.*zu_2_4.*2.0;
t128 = t63-t65;
t129 = t63-t67;
t130 = t128.^2;
t131 = t129.^2;
t132 = Wuwb.*t115.*t130.*zu_2_3.*(1.0./2.0);
t133 = Wuwb.*t110.*t117.*t130.*zu_2_3.*(1.0./2.0);
t134 = Wuwb.*t120.*t131.*zu_2_4.*(1.0./2.0);
t135 = Wuwb.*t113.*t122.*t131.*zu_2_4.*(1.0./2.0);
t136 = t56-t85-t86;
t137 = t79-t81;
t138 = t79-t83;
t139 = t137.^2;
t140 = t138.^2;
t141 = Wuwb.*t115.*t139.*zu_2_3.*(1.0./2.0);
t142 = Wuwb.*t110.*t117.*t139.*zu_2_3.*(1.0./2.0);
t143 = Wuwb.*t120.*t140.*zu_2_4.*(1.0./2.0);
t144 = Wuwb.*t113.*t122.*t140.*zu_2_4.*(1.0./2.0);
t154 = Wuwb.*t45.*t46.*zu_1_3.*(1.0./2.0);
t145 = -t48+t57-t154;
t146 = -t116-t118+t126;
t147 = Chi31-Chi41;
t148 = Chi32-Chi42;
t149 = Chi33-Chi43;
t150 = t147.^2;
t151 = t148.^2;
t152 = t149.^2;
t153 = t150+t151+t152;
t155 = t31-t35;
t156 = t155.^2;
t157 = sqrt(t153);
t158 = d_3_4-t157;
t159 = 1.0./sqrt(t153);
t160 = 1.0./t153;
t161 = Wuwb.*t156.*t160.*zu_3_4.*(1.0./2.0);
t162 = 1.0./t153.^(3.0./2.0);
t163 = Wuwb.*t156.*t158.*t162.*zu_3_4.*(1.0./2.0);
t168 = Wuwb.*t46.*t71.*zu_1_3.*(1.0./2.0);
t164 = t57-t72-t168;
t165 = t126-t132-t133;
t166 = Wgps.*zg_3.*2.0;
t167 = Wuwb.*t158.*t159.*zu_3_4.*2.0;
t169 = t65-t67;
t170 = t169.^2;
t171 = Wuwb.*t160.*t170.*zu_3_4.*(1.0./2.0);
t172 = Wuwb.*t158.*t162.*t170.*zu_3_4.*(1.0./2.0);
t175 = Wuwb.*t46.*t87.*zu_1_3.*(1.0./2.0);
t173 = t57-t88-t175;
t174 = t126-t141-t142;
t176 = t81-t83;
t177 = t176.^2;
t178 = Wuwb.*t160.*t177.*zu_3_4.*(1.0./2.0);
t179 = Wuwb.*t158.*t162.*t177.*zu_3_4.*(1.0./2.0);
t183 = Wuwb.*t50.*t51.*zu_1_4.*(1.0./2.0);
t180 = -t53+t58-t183;
t181 = -t121-t123+t127;
t182 = -t161-t163+t167;
t188 = Wuwb.*t51.*t73.*zu_1_4.*(1.0./2.0);
t184 = t58-t74-t188;
t185 = t127-t134-t135;
t186 = t167-t171-t172;
t187 = Wgps.*zg_4.*2.0;
t192 = Wuwb.*t51.*t89.*zu_1_4.*(1.0./2.0);
t189 = t58-t90-t192;
t190 = t127-t143-t144;
t191 = t167-t178-t179;
HJ = reshape([t41+t43+t48+t53+t54+t55-Wuwb.*t30.*t39.*zu_1_2.*2.0-Wuwb.*t34.*t44.*zu_1_3.*2.0-Wuwb.*t38.*t49.*zu_1_4.*2.0+Wuwb.*t25.^2.*t46.*zu_1_3.*(1.0./2.0)+Wuwb.*t26.^2.*t51.*zu_1_4.*(1.0./2.0),0.0,0.0,t91,0.0,0.0,t145,0.0,0.0,t180,0.0,0.0,0.0,t54+t55-t56-t57-t58+t69+t70+t72+t74+Wuwb.*t46.*t61.^2.*zu_1_3.*(1.0./2.0)+Wuwb.*t51.*t62.^2.*zu_1_4.*(1.0./2.0),0.0,0.0,t124,0.0,0.0,t164,0.0,0.0,t184,0.0,0.0,0.0,t54+t55-t56-t57-t58+t85+t86+t88+t90+Wuwb.*t46.*t77.^2.*zu_1_3.*(1.0./2.0)+Wuwb.*t51.*t78.^2.*zu_1_4.*(1.0./2.0),0.0,0.0,t136,0.0,0.0,t173,0.0,0.0,t189,t91,0.0,0.0,t41+t43+t54-t56+t116+t118+t121+t123+t125-Wuwb.*t110.*t114.*zu_2_3.*2.0-Wuwb.*t113.*t119.*zu_2_4.*2.0,0.0,0.0,t146,0.0,0.0,t181,0.0,0.0,0.0,t124,0.0,0.0,t54-t56+t69+t70+t125-t126-t127+t132+t133+t134+t135,0.0,0.0,t165,0.0,0.0,t185,0.0,0.0,0.0,t136,0.0,0.0,t54-t56+t85+t86+t125-t126-t127+t141+t142+t143+t144,0.0,0.0,t174,0.0,0.0,t190,t145,0.0,0.0,t146,0.0,0.0,t48+t54-t57+t116+t118-t126+t154+t161+t163+t166-Wuwb.*t158.*t159.*zu_3_4.*2.0,0.0,0.0,t182,0.0,0.0,0.0,t164,0.0,0.0,t165,0.0,0.0,t54-t57+t72-t126+t132+t133+t166-t167+t168+t171+t172,0.0,0.0,t186,0.0,0.0,0.0,t173,0.0,0.0,t174,0.0,0.0,t54-t57+t88-t126+t141+t142+t166-t167+t175+t178+t179,0.0,0.0,t191,t180,0.0,0.0,t181,0.0,0.0,t182,0.0,0.0,t53+t54-t58+t121+t123-t127+t161+t163-t167+t183+t187,0.0,0.0,0.0,t184,0.0,0.0,t185,0.0,0.0,t186,0.0,0.0,t54-t58+t74-t127+t134+t135-t167+t171+t172+t187+t188,0.0,0.0,0.0,t189,0.0,0.0,t190,0.0,0.0,t191,0.0,0.0,t54-t58+t90-t127+t143+t144-t167+t178+t179+t187+t192],[12,12]);
