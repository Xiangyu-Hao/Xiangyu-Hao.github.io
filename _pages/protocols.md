---
layout: archive
title: "Protocols"
permalink: /protocols/
author_profile: true
---
<script src="https://cdn.jsdelivr.net/npm/js-sha256@0.9.0/src/sha256.min.js"></script>
<script>
    function checkPassword() {
        // Correct SHA-256 hash for "arronhao1127"
        const correctHash = '936c53e83921ecde0dd03332463e77c8e31003f56e9c9c1598179cf0213197ef';
        
        const userPassword = prompt("请输入访问密码：");
        if (userPassword === null || userPassword === "") {
            return;
        }
        const userHash = sha256(userPassword);
        
        if (userHash === correctHash) {
            window.location.href = "../protocols/Phylogenomics.html";
        } else {
            alert("密码错误，请重试！");
        }
    }
</script>
# <a onclick="checkPassword()">Phylogenomics</a>
# <a href="../protocols/transcriptome.html" target="_blank">Transcriptome</a>
# <a href="../protocols/Molecular Convergence (conv_cal).html" target="_blank">Molecular Convergence (conv_cal)</a>
# <a href="../protocols/Relaxed Selection Test.html" target="_blank">Relaxed Selection Test</a>
